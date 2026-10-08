(* ============================================================================
   gen_grid.wl  --  generate the SHARED grid over the physical triangle
   0<x<1, 0<y<1-x, used by every dataset in this comparison pipeline (so
   c1only/coverage and weight-1/weight-2 are all compared at the exact same
   points).

   Same triangular nested-grid construction as ../no_blow_up_approach/pipeline/compare_weight.wl
   (eps-inset, per-row shrinking y-range). <N> is a TARGET, not an exact
   count: picks whichever gridN's natural triangular lattice size
   gridN(gridN+1)/2 lands closest to <N>, then keeps the FULL lattice at
   that density -- no RandomSample, no holes, full coverage of the triangle.
   Actual point count is printed and will differ slightly from <N>.

   USAGE :  wolframscript -file gen_grid.wl  <N>  <outfile>
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] < 2,
   Print["USAGE: wolframscript -file gen_grid.wl <N> <outfile>"]; Exit[1]];
targetN = ToExpression[cliArgs[[1]]];
outfile = cliArgs[[2]];

(* gridN whose triangular count gridN(gridN+1)/2 is CLOSEST to targetN
   (not just the smallest that's >= it) -- search a window around the
   analytic estimate and pick the nearest by |count-targetN|. *)
gridN0 = gridN /. Solve[gridN (gridN + 1)/2 == targetN, gridN][[-1]] // N // Round;
candidates = Select[Range[gridN0 - 3, gridN0 + 3], # >= 1 &];
gridN = candidates[[First@Ordering[Abs[#(# + 1)/2 - targetN] & /@ candidates, 1]]];

eps = 1/(2 gridN);
pts = Select[
   Flatten[Table[{xv, yv, 1 - xv - yv},
      {xv, eps, 1 - eps, (1 - 2 eps)/gridN}, {yv, eps, 1 - xv - eps, (1 - 2 eps)/gridN}], 1],
   (#[[3]] > eps) &];
Print["[grid] full triangular lattice (gridN=", gridN, ", target was ", targetN, "): ", Length[pts], " points -- entire triangle covered, no subsampling"];

Put[pts, outfile];
Print["[grid] wrote ", Length[pts], " points -> ", outfile];
