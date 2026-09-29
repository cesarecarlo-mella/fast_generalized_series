(* ============================================================================
   03_analyze.wl  --  STEP 3 (ANALYZE) : turn the raw per-component
   differences (from 01_eval_grid.wl / 02_merge_raw.wl) into the per-point
   digit-count deliverable. Cheap, single pass over already-computed data --
   re-run this as many times as you like (different chop threshold, etc.)
   WITHOUT ever re-running the expensive grid evaluation.
   ----------------------------------------------------------------------------
   Per point, per variable (native / Bernoulli):
     1. digits(component i) = relative-difference digit count, same formula
        used throughout this project (fuzzy-zero Chop[] fix included).
     2. CHOP: exclude component i from this point's aggregation if
        |ref[i]| < chopThreshold -- a component whose reference value is
        this close to zero is either exactly on, or immediately next to, a
        genuine zero-crossing of that component (confirmed directly on real
        data: e.g. component 189 of eps_2 crosses zero near x~0.26,~0.315),
        where the RELATIVE error is not a meaningful measure of convergence
        even though the series itself is fine. Excluding it here -- rather
        than letting it silently dominate the point's MIN -- is what
        produced the sharp, physically-meaningless line artifacts in the
        old single-pass pipeline.
     3. Aggregate over whatever components SURVIVE the chop:
          MIN    -- worst remaining component (the trustworthy-digit floor)
          MEDIAN -- typical remaining component (replaces the old MEAN --
                    more robust to any one remaining outlier component,
                    same reasoning as the chop itself, just one line further)

   Row (SAME shape as the old compare_*.m, just avg -> median):
     {x, y, z, corner, digN_min, digN_median, digB_min, digB_median}

       reads :  <datadir>/rawdata_<mode>_w<w>_ordp<ORDp>_ordb<ORDb>.m
       writes:  <outdir>/compare_<mode>_w<w>_ordp<ORDp>_ordb<ORDb>.m

   USAGE :  <kernel> 03_analyze.wl  <mode> <w> <ORDp> <ORDb> <datadir> [outdir] [chopThreshold]
            chopThreshold defaults to 10^-4 (absolute, on |ref[i]|)
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] < 5,
   Print["USAGE: 03_analyze.wl <mode> <w> <ORDp> <ORDb> <datadir> [outdir] [chopThreshold]"]; Exit[1]];
mode = cliArgs[[1]]; wTarget = ToExpression[cliArgs[[2]]];
ORDp = ToExpression[cliArgs[[3]]]; ORDb = ToExpression[cliArgs[[4]]];
datadir = cliArgs[[5]];
outDir = If[Length[cliArgs] >= 6, cliArgs[[6]], "."];
chopThreshold = If[Length[cliArgs] >= 7, ToExpression[cliArgs[[7]]], 10^-4];

tag = mode <> "_w" <> ToString[wTarget];
outTag = tag <> "_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDb];
inFile = FileNameJoin[{datadir, "rawdata_" <> outTag <> ".m"}];
If[! FileExistsQ[inFile], Print["[ANALYZE] ERROR: missing ", inFile]; Exit[1]];
raw = Import[inFile];
rows = raw["data"];
Print["[ANALYZE] ", tag, "  ", Length[rows], " points  chopThreshold=", chopThreshold];

ok[q_] := NumberQ[q] && Element[q, Reals];
digitOf[d_, r_] := Module[{dd, rel, res},
   If[! ok[d] || ! ok[r], Return[-2.]];
   dd = Chop[Abs[d], 10^-16];
   If[dd == 0, Return[16.]];
   rel = dd/Abs[r];
   res = Min[16., -Log10[rel + 10^-30]];
   If[ok[res], res, -2.]
];

t0 = AbsoluteTime[]; nChoppedTotal = 0; nCompTotal = 0;
(* IMPORTANT: filter to `keep` BEFORE calling digitOf, not after -- digitOf
   divides by Abs[ref], and ~136/371 weight-2 components are EXACTLY zero
   everywhere (confirmed structural fact, see 01_eval_grid.wl notes). Calling
   digitOf on those first and discarding the result afterward still means
   computing dd/Abs[0] for every one of them, at every point -- a literal
   division by zero, flooding the log with Power::infy warnings and (at
   real cluster point counts) apparently enough of them to choke the job
   entirely (confirmed: c3only/coverage order-20 jobs never finished on a
   real 400-point cluster run, while the lighter c1only/c2only modes did).
   Fix: index diffN/diffB/ref down to `keep` FIRST, so digitOf only ever
   sees components whose reference already passed the chop threshold. *)
analyzePoint[row_] := Module[
   {xv, yv, zv, c, diffN, diffB, ref, keep, diffNk0, diffBk0, refk, digNk, digBk},
   {xv, yv, zv, c, diffN, diffB, ref} = row;
   keep = Position[ref, r_ /; Abs[r] >= chopThreshold, {1}, Heads -> False][[All, 1]];
   nCompTotal += Length[ref]; nChoppedTotal += Length[ref] - Length[keep];
   If[keep === {}, keep = Range[Length[ref]]];  (* pathological fallback: nothing survived, don't crash *)
   refk = ref[[keep]]; diffNk0 = diffN[[keep]]; diffBk0 = diffB[[keep]];
   digNk = MapThread[digitOf, {diffNk0, refk}];
   digBk = MapThread[digitOf, {diffBk0, refk}];
   {xv, yv, zv, c, Min[digNk], Median[digNk], Min[digBk], Median[digBk]}
];

outRows = analyzePoint /@ rows;
Print["[ANALYZE] chopped ", nChoppedTotal, "/", nCompTotal, " (component,point) pairs  (",
      Round[100. nChoppedTotal/nCompTotal, 0.01], "%)   ", Round[AbsoluteTime[] - t0, .1], " s"];

Quiet@CreateDirectory[outDir];
outFile = FileNameJoin[{outDir, "compare_" <> outTag <> ".m"}];
Put[<|"mode" -> mode, "weight" -> wTarget, "ordp" -> ORDp, "ordb" -> ORDb, "chopThreshold" -> chopThreshold,
      "header" -> {"x", "y", "z", "corner", "digits_native_min", "digits_native_median",
                   "digits_bernoulli_min", "digits_bernoulli_median"},
      "data" -> outRows|>, outFile];
Print["[ANALYZE] wrote ", outFile, "   DONE"];
