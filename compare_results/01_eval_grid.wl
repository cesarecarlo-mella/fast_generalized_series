(* ============================================================================
   01_eval_grid.wl  --  STEP 1 (COMPUTE) : evaluate a slice of grid points,
   save the RAW per-component differences and the RAW reference vector --
   nothing about "digits", MIN/MEDIAN, or chopping small numbers happens
   here. That is all decided later, cheaply, by 03_analyze.wl reading this
   output -- so tuning the analysis (a chop threshold, a different
   aggregation statistic, ...) never requires re-running this expensive
   step again.
   ----------------------------------------------------------------------------
   Runs in the CURRENT WORKING DIRECTORY.  N processes (slice p: Mod[i-1,N]==p
   over the shared grid's point INDEX, not over components).

   For each grid point (x,y,z) in this slice:
     * mode "c<N>only" (N=1,2,3) : ALWAYS use corner N's own series,
                         evaluated everywhere in the triangle regardless of
                         where the point sits (shows corner N's own
                         convergence radius).
     * mode "coverage" : pick the NEAREST corner via the same distance-switch
                         `pickC` as ../no_blow_up_approach/pipeline/compare_weight.wl.

   Evaluate ALL 371 components, native physical variables AND Bernoulli/
   conformal variables, against a REFERENCE vector:
     * weight 1, 2 : the reference is the true exact solution (eps_1.m/
                      eps_2.m -> Log/PolyLog[2] rewrite).
     * weight 6    : no exact solution exists -- the reference is the
                      HIGHER of the two requested orders (ORDb), the lower
                      order (ORDp) is what's compared against it. Same
                      "diff vs reference" shape either way, so 03_analyze.wl
                      doesn't need to know which case it's in.

   Row: {x, y, z, corner, diffN(371), diffB(371), ref(371)}
     diffN = native_approx - ref        (SIGNED, not absolute -- keep the
     diffB = bernoulli_approx - ref      sign, 03_analyze.wl decides what to
                                          do with it)
     ref   = the reference vector itself (needed downstream to turn a
             difference into a relative measure, AND to chop components
             whose reference is too close to zero to give a meaningful
             relative comparison there)

       reads :  <griddir>/<gridfile>                 shared {x,y,z} list
                <datadir>/phys_c<c>_w<w>_ord<ORDp>.m
                <datadir>/ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>.m
                (weight 6 also reads the _ord<ORDb> / _ordp<ORDb>_ordb<ORDb> pair)
       writes:  raw_<mode>_w<w>_p<p>_<N>.m            { {i, row}, ... }
                progEVAL_<mode>_w<w>_p<p>.txt

   USAGE :  <kernel> 01_eval_grid.wl  <mode> <w> <p> <N> <gridfile> <datadir> [ORDp] [ORDb]
            mode : c1only | c2only | c3only | coverage
            weight 1,2 : ORDp = ORDb = the single order to evaluate at
            weight 6   : ORDp = lower order, ORDb = higher order (reference)
   ----------------------------------------------------------------------------
   WEIGHT 6 PERFORMANCE NOTE: weight-6 components are enormous (measured on
   real H corner-1 Bernoulli data: LeafCount from 16 up to ~119k, MEDIAN
   ~101k -- most of the 371 components, not just a couple outliers). Naive
   per-point substitution (Re@N[expr /. sub, PREC]) on an expression that
   size measured ~120 ms/point; with 371 components x up to 4 quantities
   (physLow/physHigh/berLow/berHigh) x however many grid points, that is
   hours-to-days for any real grid.

   So weight 6 (selfConvMode) uses a different strategy entirely: batch PER
   COMPONENT across this worker's whole point slice, instead of per POINT
   across all 371 components. For each component, build one tensor-
   contraction evaluator (turn the expression into a sparse monomial list
   via CoefficientRules once) and evaluate the ENTIRE slice in one dense
   matrix . vector contraction -- measured ~200x faster including the
   one-time build cost, ~75000x faster on marginal per-point cost, on the
   same real data, cross-checked against naive substitution to 1e-13 (see
   analytic.continuation/transport.boundaries/ benchmark notes -- same
   trick, first tried and verified there).

   Fully generic: auto-detects which variables carry a half-integer power
   (the sqrt/ArcTan l7/l8 sector routinely produces Sqrt[sx]-type terms
   after Bernoulli resummation) and handles exactly those, nothing hardcoded.
   Native truncations still carry literal Log[x]/Log[y]/Log[z] (never
   abstracted upstream), so those are abstracted to LGx/LGy/LGz here first to
   make them tensor-evaluable. Bernoulli files need NO such step -- they
   already use LogX/LogY/LogZ as their own native placeholder symbols (set by
   05_bernoulli.wl/06_merge_bernoulli.wl upstream) -- using the wrong name
   here once (LGx instead of LogX) silently produced a huge left-over
   symbolic mess instead of a number; caught by the spot-check below, not by
   any error message, which is exactly why that check exists.

   This changes ONLY the internal weight-6 computation strategy -- the CLI,
   the output file format/contents, and weight 1/2 (small enough that naive
   substitution is already fine) are all completely unchanged, so
   02_merge_raw.wl / 03_analyze.wl / run_compute.py need no changes at all.
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] < 6,
   Print["USAGE: 01_eval_grid.wl <mode> <w> <p> <N> <gridfile> <datadir> [ORDp] [ORDb]"]; Exit[1]];
mode    = cliArgs[[1]];
wTarget = ToExpression[cliArgs[[2]]];
part    = ToExpression[cliArgs[[3]]];
nParts  = ToExpression[cliArgs[[4]]];
gridfile = cliArgs[[5]];
datadir  = cliArgs[[6]];
ORDp = If[Length[cliArgs] >= 7, ToExpression[cliArgs[[7]]], 15];
ORDb = If[Length[cliArgs] >= 8, ToExpression[cliArgs[[8]]], 15];
fixedCorner = If[StringMatchQ[mode, "c" ~~ DigitCharacter .. ~~ "only"],
   ToExpression[StringTake[mode, {2, StringLength[mode] - 4}]], None];
If[! (MemberQ[{1, 2, 3}, fixedCorner] || mode === "coverage"),
   Print["ERROR: mode must be c1only, c2only, c3only, or coverage, got ", mode]; Exit[1]];
If[! MemberQ[{1, 2, 6}, wTarget],
   Print["ERROR: only weight 1, 2 (exact reference) or 6 (self-convergence) are supported, got ", wTarget]; Exit[1]];
selfConvMode = (wTarget === 6);

tag = mode <> "_w" <> ToString[wTarget];
cwp = tag <> "_p" <> ToString[part] <> "_" <> ToString[nParts];
progFile = "progEVAL_" <> tag <> "_p" <> ToString[part] <> ".txt";
Print["[EVAL ", cwp, "]  ORDp ", ORDp, " ORDb ", ORDb, "  datadir ", datadir, "  dir ", Directory[]];

(* ---- which corners' data this mode/slice needs -------------------------- *)
neededCorners = If[fixedCorner =!= None, {fixedCorner}, {1, 2, 3}];
physFile[c_, ord_] := FileNameJoin[{datadir, "phys_c" <> ToString[c] <> "_w" <> ToString[wTarget] <>
   "_ord" <> ToString[ord] <> ".m"}];
berFile[c_, ord_] := FileNameJoin[{datadir, "ber_c" <> ToString[c] <> "_w" <> ToString[wTarget] <>
   "_ordp" <> ToString[ord] <> "_ordb" <> ToString[ord] <> ".m"}];
loadPhysBer[c_] := If[selfConvMode,
   <|"physLow" -> Get[physFile[c, ORDp]], "physHigh" -> Get[physFile[c, ORDb]],
     "berLow" -> Get[berFile[c, ORDp]], "berHigh" -> Get[berFile[c, ORDb]]|>,
   <|"phys" -> Get[physFile[c, ORDp]],
     "ber"  -> Get[FileNameJoin[{datadir, "ber_c" <> ToString[c] <> "_w" <> ToString[wTarget] <>
                                  "_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDb] <> ".m"}]]|>];
data = Association@Table[c -> loadPhysBer[c], {c, neededCorners}];
Print["[EVAL ", cwp, "]  loaded corners ", neededCorners, "  (",
      Length[data[neededCorners[[1]]][If[selfConvMode, "physLow", "phys"]]], " components)"];

(* ---- exact reference (371-vector in x,y) -- weight 1/2 only ------------- *)
If[! selfConvMode,
 here = DirectoryName[AbsoluteFileName[
    Which[$ScriptCommandLine =!= {}, First[$ScriptCommandLine], $InputFileName =!= "", $InputFileName, True, "."]]];
 exactDir = FileNameJoin[{here, "exact_solutions"}];
 rules1 = {G[0, x] :> Log[x], G[0, y] :> Log[y], G[a_, b_] :> Log[1 - b/a]};
 rules2 = {
    G[0, x] :> Log[x], G[0, y] :> Log[y], G[0, 0, x] -> Log[x]^2/2, G[0, 0, y] -> Log[y]^2/2,
    G[0, 1, y] -> -PolyLog[2, y],
    G[1, 0, x] -> Log[1 - x] Log[x] + PolyLog[2, x],
    G[1, 0, y] -> Log[1 - y] Log[y] + PolyLog[2, y],
    G[1, 1, y] -> (1/2) Log[1 - y]^2,
    G[0, 1 - y, x] -> -PolyLog[2, -x/(-1 + y)],
    G[1 - y, 1 - y, x] -> Log[1 - x/(1 - y)]^2/2,
    G[1 - y, 0, x] -> Log[x] (Log[1 - x - y] - Log[1 - y]) + PolyLog[2, x/(1 - y)],
    G[-y, 1 - y, x] -> (Log[1 - x - y] - Log[1 - y]) Log[x + y] + PolyLog[2, 1 - x - y] - PolyLog[2, 1 - y],
    G[a_, b_] :> Log[1 - b/a]};
 exactFile = FileNameJoin[{exactDir, If[wTarget == 1, "eps_1.m", "eps_2.m"]}];
 If[! FileExistsQ[exactFile],
    Print["[EVAL ", cwp, "]  ERROR: missing ", exactFile]; Exit[1]];
 exVec = Import[exactFile] /. If[wTarget == 1, rules1, rules2];
 Print["[EVAL ", cwp, "]  exact reference loaded (", Length[exVec], " components)"],
 Print["[EVAL ", cwp, "]  self-convergence mode: N=", ORDp, " vs N=", ORDb, " (higher order = reference)"]
];

(* ---- grid slice ---------------------------------------------------------- *)
grid = Get[gridfile];
ncomp = Length[grid];
mine = Select[Range[ncomp], Mod[# - 1, nParts] == part &];
Print["[EVAL ", cwp, "]  ", Length[mine], "/", ncomp, " grid points in this slice"];

pickC[xx_, yy_, zz_] := Switch[First@Ordering[{xx, yy, zz}, -1], 1, 2, 2, 3, 3, 1];
berSub[xv_, yv_, zv_] := {LogX -> Log[xv], LogY -> Log[yv], LogZ -> Log[zv],
   sx -> -Log[1 - xv], sy -> -Log[1 - yv], sz -> -Log[1 - zv], x -> xv, y -> yv, z -> zv};
PREC = 20;  (* plenty of headroom above the 16-digit cap any downstream metric will use *)
okV[v_] := VectorQ[v, NumberQ[#] && Element[#, Reals] &];

(* row: {x,y,z,corner, diffN(371), diffB(371), ref(371)} -- SIGNED diffs *)
evalPoint[xv_, yv_, zv_] := Module[{c, sub, phys, ber, ref, appN, appB, diffN, diffB},
   c = If[fixedCorner =!= None, fixedCorner, pickC[xv, yv, zv]];
   sub = {x -> xv, y -> yv, z -> zv};
   phys = data[c]["phys"]; ber = data[c]["ber"];
   ref = Re@N[exVec /. sub, PREC];
   appN = Re@N[phys /. sub, PREC];
   appB = Re@N[ber /. berSub[xv, yv, zv], PREC];
   diffN = appN - ref; diffB = appB - ref;
   If[! (okV[diffN] && okV[diffB] && okV[ref]), Return[$Failed]];
   {N[xv], N[yv], N[zv], c, diffN, diffB, ref}
];

(* ============================================================================
   FAST PATH (weight 6 / selfConvMode only) : tensor-contraction batch
   evaluator, per component, over a whole slice of grid points at once.
   ============================================================================ *)

buildTensorEvaluator[expr_, vars_List] := Module[
   {halfVars, subVars, sqrtSubs, expr2, rules, exps, coeffs, maxE},
  halfVars = Select[vars, ! FreeQ[expr, Power[#, e_] /; ! IntegerQ[e]] &];
  subVars = Table[If[MemberQ[halfVars, v], Unique["q"], v], {v, vars}];
  sqrtSubs = Thread[halfVars -> (Extract[subVars, Position[vars, #]][[1]]^2 & /@ halfVars)];
  expr2 = PowerExpand[Expand[expr /. sqrtSubs],
     Assumptions -> Thread[subVars[[Flatten[Position[vars, #] & /@ halfVars]]] >= 0]];
  rules = CoefficientRules[expr2, subVars];
  exps = rules[[All, 1]];
  coeffs = N[rules[[All, 2]]];
  maxE = If[exps === {}, ConstantArray[0, Length[vars]], Max /@ Transpose[exps]];
  <|"vars" -> vars, "halfVars" -> halfVars, "exps" -> exps, "coeffs" -> coeffs, "maxE" -> maxE|>
];

(* valTable: {ngrid x nvars} numeric matrix, values for `vars` IN ORDER, in
   the ORIGINAL (non-sqrt-substituted) variables -- Sqrt taken internally for
   whichever of `vars` were flagged half-integer. *)
evalTensorBatch[spec_Association, valTable_] := Module[
   {vars = spec["vars"], halfVars = spec["halfVars"], maxE = spec["maxE"],
    exps = spec["exps"], coeffs = spec["coeffs"], ngrid = Length[valTable], pows, monMat},
  If[exps === {}, Return[ConstantArray[0., ngrid]]];
  pows = Table[
     Module[{colVals = valTable[[All, i]]},
        If[MemberQ[halfVars, vars[[i]]], colVals = Sqrt[colVals]];
        Outer[Power, colVals, Range[0, maxE[[i]]]]],
     {i, Length[vars]}];
  monMat = Transpose[Table[
     Times @@ Table[pows[[i]][[All, e[[i]] + 1]], {i, Length[vars]}],
     {e, exps}]];
  monMat . coeffs
];

physVars = {x, y, z, LGx, LGy, LGz};
(* Bernoulli output files already use LogX/LogY/LogZ as their OWN native
   placeholder symbols (set by 05_bernoulli.wl/06_merge_bernoulli.wl in the
   blow_up_form_2 pipeline) -- no abstraction needed/applied here, unlike the
   native side above. Confirmed directly: `vars in biggest [component]` of a
   real ber_c1_w6_...m file returned {LogX, LogY, Pi, sx, sy}, never LGx/LGy. *)
berVars  = {sx, sy, sz, LogX, LogY, LogZ};

(* one-point independent naive cross-check, printed once per corner, so a
   silent bug in the tensor machinery can never pass unnoticed *)
spotCheckCorner[c_, xv0_, yv0_, zv0_] := Module[
   {pL = data[c]["physLow"], bL = data[c]["berLow"], k, fastN, naiveN, fastB, naiveB, specN, specB},
  k = Max[1, Round[Length[pL]/2]];  (* a "typical", not-necessarily-trivial component *)
  specN = buildTensorEvaluator[pL[[k]] /. {Log[x] -> LGx, Log[y] -> LGy, Log[z] -> LGz}, physVars];
  fastN = First[evalTensorBatch[specN, {{xv0, yv0, zv0, Log[xv0], Log[yv0], Log[zv0]}}]];
  naiveN = Re@N[pL[[k]] /. {x -> xv0, y -> yv0, z -> zv0}, 20];
  specB = buildTensorEvaluator[bL[[k]], berVars];
  fastB = First[evalTensorBatch[specB, {{-Log[1 - xv0], -Log[1 - yv0], -Log[1 - zv0], Log[xv0], Log[yv0], Log[zv0]}}]];
  naiveB = Re@N[bL[[k]] /. berSub[xv0, yv0, zv0], 20];
  If[Abs[fastN - naiveN] > 10^-8 Max[1, Abs[naiveN]] || Abs[fastB - naiveB] > 10^-8 Max[1, Abs[naiveB]],
     Print["[EVAL ", cwp, "]  !! FAST-PATH SPOT-CHECK FAILED  corner ", c, " component ", k,
           "  native fast=", fastN, " naive=", naiveN, "  bernoulli fast=", fastB, " naive=", naiveB],
     Print["[EVAL ", cwp, "]  [check] corner ", c, " component ", k,
           " fast-path verified against naive substitution (native diff ",
           ScientificForm[N[Abs[fastN - naiveN]], 2], ", bernoulli diff ",
           ScientificForm[N[Abs[fastB - naiveB]], 2], ")"]];
];

(* Optional precomputed cache, from 00_build_coeffcache.wl -- see that file's
   header. If present, skips the Expand+PowerExpand+CoefficientRules step
   entirely (the dominant cost, ~16 min/corner) for THIS (corner,ORDp,ORDb),
   loading its {exps,coeffs} specs directly instead. Purely additive: if the
   file isn't there, fastEvalCornerSlice below builds specs fresh exactly as
   before this change -- nothing required, nothing breaks either way. Only
   applies in selfConvMode (weight 6); checked once, outside the point loop. *)
coeffCache = If[selfConvMode,
   Association@Table[
      Module[{f = FileNameJoin[{datadir, "coeffcache_c" <> ToString[c] <> "_w6_ordp" <>
                 ToString[ORDp] <> "_ordb" <> ToString[ORDb] <> ".m"}]},
         If[FileExistsQ[f],
          (
           Print["[EVAL ", cwp, "]  corner ", c, "  using coeff cache ", f];
           c -> Get[f]
          ),
           c -> None]],
      {c, neededCorners}],
   Association[]];

(* evaluate an entire corner's worth of points (idxList, indices into `grid`)
   for ALL components at once -- returns {{i,row},...} in the SAME format
   evalPointSelfConv would have produced pointwise.

   (A shared-power-table variant -- building each variable's Outer[Power,...]
   table once instead of once per component -- was tried and measured on
   real weight-6 data: only a 1.06x speedup at 2000 points/10 components,
   because the dominant cost is Expand+PowerExpand+CoefficientRules, which
   is inherently per-component and can't be shared; the power-table build it
   targeted was already a small fraction of the total. Not worth the added
   complexity/fallback-logic risk for that gain -- reverted. The coeffCache
   above targets the SAME dominant cost properly: since it depends only on
   the component's own expression, not the grid, precomputing it ONCE
   (00_build_coeffcache.wl) and loading it here removes the redundant work
   across every parallel point-slicing worker, instead of trying to share
   less-expensive parts within one worker.) *)
fastEvalCornerSlice[c_, idxList_] := Module[
   {ptsC, xv, yv, zv, pL, pH, bL, bH, ncompLocal, valTableN, valTableB,
    diffNmat, diffBmat, refmat, tcomp, rows, cached = coeffCache[c]},
  If[idxList === {}, Return[{}]];
  ptsC = grid[[idxList]];
  xv = ptsC[[All, 1]]; yv = ptsC[[All, 2]]; zv = ptsC[[All, 3]];
  valTableN = Transpose[{xv, yv, zv, Log[xv], Log[yv], Log[zv]}];
  valTableB = Transpose[{-Log[1 - xv], -Log[1 - yv], -Log[1 - zv], Log[xv], Log[yv], Log[zv]}];
  spotCheckCorner[c, xv[[1]], yv[[1]], zv[[1]]];
  pL = data[c]["physLow"]; pH = data[c]["physHigh"];
  bL = data[c]["berLow"]; bH = data[c]["berHigh"];
  ncompLocal = Length[pL];
  tcomp = AbsoluteTime[];
  {diffNmat, refmat, diffBmat} = Transpose@Table[
     Module[{eLn, eHn, specLn, specHn, vLn, vHn, specLb, specHb, vLb, vHb},
        If[cached =!= None,
         (
          specLn = cached["physLow"][[k]]; specHn = cached["physHigh"][[k]];
          specLb = cached["berLow"][[k]]; specHb = cached["berHigh"][[k]]
         ),
         (
          eLn = pL[[k]] /. {Log[x] -> LGx, Log[y] -> LGy, Log[z] -> LGz};
          eHn = pH[[k]] /. {Log[x] -> LGx, Log[y] -> LGy, Log[z] -> LGz};
          specLn = buildTensorEvaluator[eLn, physVars];
          specHn = buildTensorEvaluator[eHn, physVars];
          specLb = buildTensorEvaluator[bL[[k]], berVars];
          specHb = buildTensorEvaluator[bH[[k]], berVars];
         )];
        vLn = Re[N[evalTensorBatch[specLn, valTableN]]];
        vHn = Re[N[evalTensorBatch[specHn, valTableN]]];
        vLb = Re[N[evalTensorBatch[specLb, valTableB]]];
        vHb = Re[N[evalTensorBatch[specHb, valTableB]]];
        If[Mod[k, 25] == 0 || k == ncompLocal,
           Print["[EVAL ", cwp, "]  corner ", c, " fast-path component ", k, "/", ncompLocal,
                 "  ", Round[AbsoluteTime[] - tcomp, 0.1], " s"]];
        {vLn - vHn, vHn, vLb - vHb}
     ], {k, ncompLocal}];
  (* diffNmat/refmat/diffBmat : ncompLocal x npts -- transpose to per-point rows *)
  rows = Table[
     {N[xv[[j]]], N[yv[[j]]], N[zv[[j]]], c,
      diffNmat[[All, j]], diffBmat[[All, j]], refmat[[All, j]]},
     {j, Length[idxList]}];
  Transpose[{idxList, rows}]
];

t0 = AbsoluteTime[]; out = {}; nd = 0; nn = Length[mine];

If[selfConvMode,
 (
  (* -------------------- FAST PATH : weight 6 -------------------- *)
  Module[{byCorner},
   byCorner = If[fixedCorner =!= None,
      <|fixedCorner -> mine|>,
      GroupBy[mine, pickC @@ grid[[#]] &]];
   Do[
    If[KeyExistsQ[byCorner, c],
     (
      Print["[EVAL ", cwp, "]  fast path: corner ", c, "  ", Length[byCorner[c]], " points"];
      out = Join[out, fastEvalCornerSlice[c, byCorner[c]]];
      nd = Length[out];
      Put[ToString[nd] <> " " <> ToString[nn], progFile];
      Print["[EVAL ", cwp, "]  corner ", c, " done  ", nd, "/", nn, "  ",
            Round[AbsoluteTime[] - t0, 0.1], " s"]
     )],
    {c, {1, 2, 3}}];
  ]
 ),
 (
  (* ---------------- ORIGINAL per-point PATH : weight 1,2 (unchanged) ---------------- *)
  Do[
     tc = AbsoluteTime[];
     pt = grid[[i]];
     r = evalPoint @@ pt;
     dt = AbsoluteTime[] - tc;
     AppendTo[out, {i, r}];
     nd++;
     If[dt > 20, Print["[EVAL ", cwp, "]  slow point ", i, " : ", Round[dt, 0.1], " s"]];
     If[Mod[nd, 25] == 0 || nd == nn,
        Put[ToString[nd] <> " " <> ToString[nn], progFile];
        Print["[EVAL ", cwp, "]  ", nd, "/", nn, "  (pt ", i, ")  ",
              Round[AbsoluteTime[] - t0, 0.1], " s"];
     ],
     {i, mine}
  ]
 )];

Put[out, "raw_" <> cwp <> ".m"];
Print["[EVAL ", cwp, "] DONE  ", Round[AbsoluteTime[] - t0, .1], " s"];
