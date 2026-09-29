(* ============================================================================
   chop_phys_ber.wl  --  derive a LOWER-order phys_/ber_ pair from an
   already-computed higher-order one, by truncation alone (no recursion, no
   Bernoulli re-expansion) -- same "chop from a master" principle as
   pipeline_master/chop_letters.wl, applied to the physical-variable /
   Bernoulli deliverables instead of the raw letter expansions.

   Exact, not approximate -- plain Select on stored exponents, same
   discard-is-safe reasoning used everywhere else in this project:

     phys_ : BOTH the corner's SMALL and BIG variables (x&y, z&y, or x&z)
             get re-truncated to <= ordTarget. NOTE: this used to chop only
             BIG, leaving SMALL fixed at the original recursion's ORDt --
             matching 03_to_physical_vars.wl's own truncBIG[], which exists
             to drop rectangle-cutoff pollution, NOT to define a genuine
             two-variable truncation order. That meant an "ord18 vs ord19"
             self-convergence check (01_eval_grid.wl) was structurally BLIND
             to any non-convergence in the SMALL variable, since SMALL never
             changed between the two files being compared -- confirmed on
             real weight-6 data: several components showed <1 digit of
             genuine (x,y)-joint convergence at a point where the old
             BIG-only check reported ~13 digits. Chopping BOTH here now
             (matching what ber_ already correctly did) fixes that -- this
             is a correctness fix, not an alternate mode; "ord" now means
             the same thing for phys_ that it always meant for ber_.
     ber_  : BOTH of the corner's s-variables (sx&sy, sy&sz, or sx&sz) get
             independently re-truncated to <= ordTarget -- matches exactly
             what 05_bernoulli.wl's own per-variable Series truncation does.
             Unchanged by this fix.

   USAGE :  wolframscript -file chop_phys_ber.wl <corner> <weight> <ordMax> <ordTarget> <datadir> <outdir>
            reads  <datadir>/phys_c<c>_w<w>_ord<ordMax>.m
                   <datadir>/ber_c<c>_w<w>_ordp<ordMax>_ordb<ordMax>.m
            writes <outdir>/phys_c<c>_w<w>_ord<ordTarget>.m
                   <outdir>/ber_c<c>_w<w>_ordp<ordTarget>_ordb<ordTarget>.m
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] < 6,
   Print["USAGE: chop_phys_ber.wl <corner> <weight> <ordMax> <ordTarget> <datadir> <outdir>"]; Exit[1]];
corner = ToExpression[cliArgs[[1]]];
wTarget = ToExpression[cliArgs[[2]]];
ordMax = ToExpression[cliArgs[[3]]];
ordTarget = ToExpression[cliArgs[[4]]];
datadir = cliArgs[[5]];
outdir = cliArgs[[6]];
If[ordTarget > ordMax,
   Print["ERROR: ordTarget (", ordTarget, ") > ordMax (", ordMax, ") -- nothing to chop from"]; Exit[1]];

bigOf = <|1 -> y, 2 -> y, 3 -> z|>;
smallOf = <|1 -> x, 2 -> z, 3 -> x|>;
sVarsOf = <|1 -> {sx, sy}, 2 -> {sy, sz}, 3 -> {sx, sz}|>;
BIG = bigOf[corner];
SMALL = smallOf[corner];
sVars = sVarsOf[corner];

cutBoth[e_, v1_, v2_, ord_] := Module[{u = Expand[e]},
   Total@Select[If[Head[u] === Plus, List @@ u, {u}],
      (Exponent[#, v1] <= ord && Exponent[#, v2] <= ord) &]];

t0 = AbsoluteTime[];
physIn = FileNameJoin[{datadir, "phys_c" <> ToString[corner] <> "_w" <> ToString[wTarget] <>
   "_ord" <> ToString[ordMax] <> ".m"}];
berIn = FileNameJoin[{datadir, "ber_c" <> ToString[corner] <> "_w" <> ToString[wTarget] <>
   "_ordp" <> ToString[ordMax] <> "_ordb" <> ToString[ordMax] <> ".m"}];
If[! FileExistsQ[physIn], Print["ERROR: missing ", physIn]; Exit[1]];
If[! FileExistsQ[berIn], Print["ERROR: missing ", berIn]; Exit[1]];

phys = Get[physIn];
ber = Get[berIn];
Print["[chop] loaded phys (", Length[phys], " comps) + ber (", Length[ber], " comps)  ",
      Round[AbsoluteTime[] - t0, .1], " s"];

t1 = AbsoluteTime[];
physOut = cutBoth[#, SMALL, BIG, ordTarget] & /@ phys;
berOut = cutBoth[#, sVars[[1]], sVars[[2]], ordTarget] & /@ ber;
Print["[chop] truncated  ", Round[AbsoluteTime[] - t1, .1], " s"];

Quiet@CreateDirectory[outdir];
physOutFile = FileNameJoin[{outdir, "phys_c" <> ToString[corner] <> "_w" <> ToString[wTarget] <>
   "_ord" <> ToString[ordTarget] <> ".m"}];
berOutFile = FileNameJoin[{outdir, "ber_c" <> ToString[corner] <> "_w" <> ToString[wTarget] <>
   "_ordp" <> ToString[ordTarget] <> "_ordb" <> ToString[ordTarget] <> ".m"}];
Put[physOut, physOutFile];
Put[berOut, berOutFile];
Print["[chop] wrote ", physOutFile];
Print["[chop] wrote ", berOutFile];
Print["DONE  ", Round[AbsoluteTime[] - t0, .1], " s total"];
