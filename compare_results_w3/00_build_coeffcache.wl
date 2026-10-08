(* ============================================================================
   00_build_coeffcache.wl  --  ONE-TIME precompute for the weight-6 fast path
   in 01_eval_grid.wl.
   ----------------------------------------------------------------------------
   Turns each of the 371 components' phys/ber expressions into their sparse
   {exps,coeffs} monomial representation (buildTensorEvaluator's own output
   format -- see 01_eval_grid.wl) and saves it to disk. This is the DOMINANT
   cost of the weight-6 fast path (measured: ~16 min/corner, ~1484
   Expand+PowerExpand+CoefficientRules calls), and it depends ONLY on the
   component's own expression -- never on the grid or how many points are
   being evaluated. So it's identical every time, yet every parallel
   point-slicing worker in 01_eval_grid.wl currently rebuilds it from
   scratch. Caching it once here means N workers (any grid, any point count)
   all just Get this file instead of redoing that ~16 min of symbolic work.

   Run this ONCE per (corner, ORDp, ORDb) combination you intend to compare,
   before run_compute.py / 01_eval_grid.wl for weight 6. If the cache file
   is present, 01_eval_grid.wl uses it automatically (near-instant); if not,
   it builds fresh exactly as before this change -- nothing breaks if you
   skip this step, it's purely an optional speedup.

   USAGE :  <kernel> 00_build_coeffcache.wl  <c> <ORDp> <ORDb> <datadir> [outdir]
       c       : corner (1,2,3)
       ORDp    : lower order (weight 6 self-convergence)
       ORDb    : higher order (reference)
       datadir : dir with phys_c<c>_w6_ord<ORDp/ORDb>.m,
                          ber_c<c>_w6_ordp<ORDp/ORDb>_ordb<ORDp/ORDb>.m
       outdir  : where to write the cache file (default datadir)
   writes:  <outdir>/coeffcache_c<c>_w6_ordp<ORDp>_ordb<ORDb>.m
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];
If[Length[cliArgs] < 4,
   Print["USAGE: 00_build_coeffcache.wl <c> <ORDp> <ORDb> <datadir> [outdir]"]; Exit[1]];
corner  = ToExpression[cliArgs[[1]]];
ORDp    = ToExpression[cliArgs[[2]]];
ORDb    = ToExpression[cliArgs[[3]]];
datadir = cliArgs[[4]];
outdir  = If[Length[cliArgs] >= 5, cliArgs[[5]], datadir];
Quiet@CreateDirectory[outdir];

physVars = {x, y, z, LGx, LGy, LGz};
berVars  = {sx, sy, sz, LogX, LogY, LogZ};

(* identical to buildTensorEvaluator in 01_eval_grid.wl -- output format
   (vars/halfVars/exps/coeffs/maxE) is Put/Get-safe as-is: exps/coeffs are
   pure numbers, halfVars/vars are plain global symbols (sx, LogX, ...),
   never the Unique["q"] symbols used internally during CoefficientRules. *)
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

t0 = AbsoluteTime[];
pL = Get[FileNameJoin[{datadir, "phys_c" <> ToString[corner] <> "_w6_ord" <> ToString[ORDp] <> ".m"}]];
pH = Get[FileNameJoin[{datadir, "phys_c" <> ToString[corner] <> "_w6_ord" <> ToString[ORDb] <> ".m"}]];
bL = Get[FileNameJoin[{datadir, "ber_c" <> ToString[corner] <> "_w6_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDp] <> ".m"}]];
bH = Get[FileNameJoin[{datadir, "ber_c" <> ToString[corner] <> "_w6_ordp" <> ToString[ORDb] <> "_ordb" <> ToString[ORDb] <> ".m"}]];
ncomp = Length[pL];
Print["[cache] corner ", corner, "  ORDp ", ORDp, " ORDb ", ORDb, "  ", ncomp, " components  loaded  ",
      Round[AbsoluteTime[] - t0, .1], " s"];

t0 = AbsoluteTime[];
physLowSpecs = {}; physHighSpecs = {}; berLowSpecs = {}; berHighSpecs = {};
Do[
 Module[{eLn, eHn},
   eLn = pL[[k]] /. {Log[x] -> LGx, Log[y] -> LGy, Log[z] -> LGz};
   eHn = pH[[k]] /. {Log[x] -> LGx, Log[y] -> LGy, Log[z] -> LGz};
   AppendTo[physLowSpecs, buildTensorEvaluator[eLn, physVars]];
   AppendTo[physHighSpecs, buildTensorEvaluator[eHn, physVars]];
   AppendTo[berLowSpecs, buildTensorEvaluator[bL[[k]], berVars]];
   AppendTo[berHighSpecs, buildTensorEvaluator[bH[[k]], berVars]];
   If[Mod[k, 25] == 0 || k == ncomp,
      Print["[cache] component ", k, "/", ncomp, "  ", Round[AbsoluteTime[] - t0, 0.1], " s"]];
 ],
 {k, ncomp}];

outfile = FileNameJoin[{outdir, "coeffcache_c" <> ToString[corner] <> "_w6_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDb] <> ".m"}];
Put[<|"physLow" -> physLowSpecs, "physHigh" -> physHighSpecs, "berLow" -> berLowSpecs, "berHigh" -> berHighSpecs|>, outfile];
Print["[cache] wrote ", outfile, "  (build ", Round[AbsoluteTime[] - t0, .1], " s, ",
      Round[FileByteCount[outfile]/10.^6, 0.1], " MB)"];
