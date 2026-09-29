(* ============================================================================
   06_merge_bernoulli.wl  --  STEP 6 : merge the Bernoulli-series slices
   ----------------------------------------------------------------------------
   Runs in the CURRENT WORKING DIRECTORY.  1 process.  Trivial / fast.

       reads :  ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>_p<p>_<N>.m   for p = 0..N-1
                                                                   ({ {i, r_i} })
       writes:  ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>.m             the deliverable
                (371-vector, double series in sx/sy/sz depending on corner,
                LogX/LogY/LogZ placeholders for Log[x]/Log[y]/Log[z])

   USAGE :  <kernel> 06_merge_bernoulli.wl  <c> <w> <ORDp> <ORDb> <N> [outdir]
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] >= 5,
   {corner, wTarget, ORDp, ORDb, nParts} = ToExpression /@ Take[cliArgs, 5],
   {corner, wTarget, ORDp, ORDb, nParts} = {1, 1, 15, 15, 1}
];
outDir = If[Length[cliArgs] >= 6, cliArgs[[6]], "."];
cwpb = "c" <> ToString[corner] <> "_w" <> ToString[wTarget] <>
       "_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDb];

acc = <||>;
Do[
   f = "ber_" <> cwpb <> "_p" <> ToString[p] <> "_" <> ToString[nParts] <> ".m";
   If[! FileExistsQ[f], Print["[06M] ABORT: missing ", f]; Exit[1]];
   Do[acc[e[[1]]] = e[[2]], {e, Get[f]}],
   {p, 0, nParts - 1}
];
ber = Values[KeySort[acc]];
Print["[06M] merged ", Length[ber], " components"];
Quiet@CreateDirectory[outDir];
outFile = FileNameJoin[{outDir, "ber_" <> cwpb <> ".m"}];
Put[ber, outFile];
Print["[06M] wrote ", outFile, "   (", ByteCount[ber], " b)   DONE"];
