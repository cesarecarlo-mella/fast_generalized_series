(* ============================================================================
   02_merge_raw.wl  --  STEP 2 (COMPUTE) : merge the raw-evaluation slices
   from 01_eval_grid.wl into one deliverable per mode/weight/order(pair).
   ----------------------------------------------------------------------------
   Runs in the CURRENT WORKING DIRECTORY.  1 process.

       reads :  raw_<mode>_w<w>_p<p>_<N>.m   for p = 0..N-1   ({ {i, row}, ... })
       writes:  rawdata_<mode>_w<w>.m         the deliverable:
                <| "mode"->mode, "weight"->w, "ordp"->ORDp, "ordb"->ORDb,
                   "header"->{"x","y","z","corner","diffN","diffB","ref"},
                   "data"-> { row, ... } |>          (rows in original grid order)

   USAGE :  <kernel> 02_merge_raw.wl  <mode> <w> <ORDp> <ORDb> <N> [outdir]
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] < 5,
   Print["USAGE: 02_merge_raw.wl <mode> <w> <ORDp> <ORDb> <N> [outdir]"]; Exit[1]];
mode = cliArgs[[1]]; wTarget = ToExpression[cliArgs[[2]]];
ORDp = ToExpression[cliArgs[[3]]]; ORDb = ToExpression[cliArgs[[4]]];
nParts = ToExpression[cliArgs[[5]]];
outDir = If[Length[cliArgs] >= 6, cliArgs[[6]], "."];
tag = mode <> "_w" <> ToString[wTarget];
outTag = tag <> "_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDb];

acc = <||>;
Do[
   f = "raw_" <> tag <> "_p" <> ToString[p] <> "_" <> ToString[nParts] <> ".m";
   If[! FileExistsQ[f], Print["[02M] ABORT: missing ", f]; Exit[1]];
   Do[
      If[e[[2]] =!= $Failed, acc[e[[1]]] = e[[2]]],
      {e, Get[f]}],
   {p, 0, nParts - 1}
];
rows = Values[KeySort[acc]];
Print["[02M] merged ", Length[rows], " points"];
Quiet@CreateDirectory[outDir];
outFile = FileNameJoin[{outDir, "rawdata_" <> outTag <> ".m"}];
Put[<|"mode" -> mode, "weight" -> wTarget, "ordp" -> ORDp, "ordb" -> ORDb,
      "header" -> {"x", "y", "z", "corner", "diffN", "diffB", "ref"},
      "data" -> rows|>, outFile];
Print["[02M] wrote ", outFile, "   (", ByteCount[rows], " b)   DONE"];
