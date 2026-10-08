(* ::Package:: *)

(* Triangle convergence maps -- reads output of run_compute.py + run_analyze.py
   (the two-stage pipeline: run_compute.py does the expensive grid evaluation
   and saves raw differences, run_analyze.py cheaply turns those into the
   compare_*.m files this script reads, with a tunable --chop threshold for
   excluding near-zero-reference components from MIN/MEDIAN).

   HOW TO USE: open a NEW blank notebook in Mathematica, paste this ENTIRE
   file into ONE cell, click into it, press Shift+Enter. Plain code, no
   file-opening tricks needed.

   Reads compare_<mode>_w<w>_ordp<O>_ordb<O>.m from resultsDir (mode in
   {c1only,c2only,c3only,coverage}). Each row is
   {x,y,z,corner,digits_native_min,digits_native_median,
    digits_bernoulli_min,digits_bernoulli_median} -- MIN is the worst-converged
   of all (non-chopped) components (the trustworthy-digit estimate), MEDIAN is
   the median over all (non-chopped) components (a less pessimistic summary,
   more robust to any one remaining outlier than a plain mean would be).
   Both are plotted.

   THREE checks, each its own bare-expression output (no trailing ";"), so
   they print one after another when the cell is evaluated. "x,y,z" below
   means the plain physical (native) kinematic variables, as opposed to the
   Bernoulli/conformal (sx,sy,sz) ones -- every panel plots BOTH min and
   median, for both:
     A. weight 2 vs exact, single order (mainOrd) -- cornerMode vs coverage,
        x,y,z vs Bernoulli, MIN vs MEDIAN (one 4x2 grid of triangles).
     B. weight 2 vs exact, evolution across w2Orders -- cornerMode vs
        coverage, x,y,z vs Bernoulli, min & median (8-row grid).
     C. weight 6 self-convergence, the two adjacent-order pairs in w6Pairs
        (18-19 then 19-20 by default) -- if pair 2's digits are higher (and
        more uniform) than pair 1's, the series is converging as N grows.
        cornerMode vs coverage, x,y,z vs Bernoulli, min & median (8-row grid).

   Each check is ALSO Export'ed as a 300dpi PNG to figures/ (local disk only,
   nothing uploaded anywhere) -- ready to \includegraphics in a paper. PNG
   not PDF: these are smooth interpolated color fields (ListDensityPlot),
   and exporting that as vector PDF blows up to tens of MB from mesh
   polygons: standard practice for this kind of continuous heatmap is a
   rasterized image at print resolution, same as any journal figure. *)

SetDirectory["/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/series.all/compare_results"];
resultsDirW2 = "results";           (* <-- edit: where weight-1/2 compare_*.m / grid_*.m live *)
resultsDirW6 = "results_w6/results"; (* <-- edit: where weight-6 compare_*.m / grid_*.m live *)
figDir = "figures"; Quiet@CreateDirectory[figDir];
figDPI = 150;  (* <-- edit: lower = smaller files/coarser, e.g. 96 for quick drafts, 300 for max print quality *)
cornerMode = "c1only";     (* <-- edit: c1only | c2only | c3only, paired with coverage in every check *)
w2Orders = {12, 16, 20};   (* <-- edit: orders you ran run_compute.py/run_analyze.py at *)
mainOrd = Last[w2Orders];  (* <-- order used for check A's single-order snapshot *)
w6Pairs = {{18, 19}, {19, 20}};  (* <-- edit: CONSECUTIVE pairs produced by --w6-orders *)

(* which two physical variables go to 0 at this corner (see 01_eval_grid.wl's
   own "corner 1: IN=y OU=x  corner 2: IN=z OU=y  corner 3: IN=x OU=z"): *)
cornerLabelOf = <|"c1only" -> "(x,y)~(0,0)", "c2only" -> "(y,z)~(0,0)",
   "c3only" -> "(x,z)~(0,0)", "coverage" -> "coverage"|>;
cornerLabel = cornerLabelOf[cornerMode];

dirFor[w_] := If[w == 6, resultsDirW6, resultsDirW2];
compareFile[mode_, w_, ordp_, ordb_] := FileNameJoin[{dirFor[w], "compare_" <> mode <> "_w" <> ToString[w] <>
   "_ordp" <> ToString[ordp] <> "_ordb" <> ToString[ordb] <> ".m"}];
haveData[mode_, w_, ordp_, ordb_] := FileExistsQ[compareFile[mode, w, ordp, ordb]];
load[mode_, w_, ordp_, ordb_] := Import[compareFile[mode, w, ordp, ordb]]["data"];

(* column layout: 1=x 2=y 3=z 4=corner 5=digN_min 6=digN_median 7=digB_min 8=digB_median *)
colOf[native_, metric_] := Which[native && metric == "min", 5, native && metric == "median", 6,
   !native && metric == "min", 7, True, 8];

cf[q_] := Blend[{{-2, RGBColor[0.75, 0.1, 0.1]}, {2, RGBColor[0.75, 0.1, 0.1]},
   {3, RGBColor[0.9, 0.35, 0.05]}, {4, RGBColor[0.95, 0.6, 0.05]}, {5, RGBColor[1, 0.85, 0.15]},
   {6, RGBColor[0.75, 0.88, 0.35]}, {8, RGBColor[0.45, 0.75, 0.35]}, {10, RGBColor[0.2, 0.6, 0.25]},
   {12, RGBColor[0.1, 0.45, 0.15]}, {16, RGBColor[0.02, 0.25, 0.05]}}, Clip[q, {-2, 16}]];
legend = BarLegend[{cf, {-2, 16}}, {-2, 2, 3, 4, 5, 6, 8, 10, 12, 16}, LegendLabel -> "digits"];
withLegend[grid_] := Row[{legend, Spacer[15], grid}];  (* legend to the LEFT of the grid, every export *)
(* Smooth-looking field WITHOUT Mathematica's scattered-data triangulation
   (ListDensityPlot/Interpolation): gen_grid.wl's points sit on a near-exact
   REGULAR triangular lattice, and that regularity makes Mathematica's
   Delaunay/FEM mesh degenerate (confirmed: raising MaxPlotPoints made it
   WORSE -- more torn/patchy, not smoother; jittering the points didn't fix
   it either). Since we KNOW the true lattice structure, draw it directly:
   one filled Rectangle per data point, sized to exactly the grid step (tiny
   1.02x overlap so there's no antialiasing seam between neighbors) -- no
   interpolation, no mesh, so it can't degenerate. ColorFunctionScaling has
   no meaning here since cf[] is called directly on the RAW digit value. *)
gridStep[d_] := Min[Differences[Sort[DeleteDuplicates[Round[d[[All, 1]], 10^-6]]]]];
(* lbl (top, via PlotLabel): expansion corner + min/median ONLY.
   insetText (in-plot, via a framed box sitting in the empty white triangle
   at top-right, x+y>1): which representation (x,y,z / Bernoulli) + the
   max_power info -- out of the title, inside the panel instead. *)
mkTri[d_, col_, lbl_, insetText_ : ""] := Module[{hw = gridStep[d]/2 * 1.02},
   Graphics[Prepend[
      Map[{cf[#[[col]]], Rectangle[{#[[1]], #[[2]]} - hw, {#[[1]], #[[2]]} + hw]} &, d],
      EdgeForm[None]],
   Frame -> True, PlotLabel -> lbl, AspectRatio -> 1, ImageSize -> 220,
   PlotRange -> {{0, 1}, {0, 1}}, FrameLabel -> {"x", "y"},
   Epilog -> {GrayLevel[0.45], Line[{{0, 0}, {1, 0}, {0, 1}, {0, 0}}],
      If[insetText =!= "",
         Inset[Framed[Style[insetText, 9], Background -> White,
            FrameStyle -> GrayLevel[0.7], FrameMargins -> 3], {0.72, 0.72}],
         {}]}]];

(* ---- check A: weight 2 vs exact, single order mainOrd -- cornerMode & coverage, x,y,z & Bernoulli, MIN & MEDIAN -- *)
dCornA = load[cornerMode, 2, mainOrd, mainOrd];
dCovA = load["coverage", 2, mainOrd, mainOrd];
titleCornA = cornerLabel <> ", max_pow=" <> ToString[mainOrd];
titleCovA = "coverage, max_pow=" <> ToString[mainOrd];
checkA = Column[{Style["weight 2  max_power=" <> ToString[mainOrd] <> "  vs exact  --  " <> cornerLabel <>
    " & coverage, x,y,z & Bernoulli, min & median", 14, Bold],
   withLegend[Grid[{
     {mkTri[dCornA, 5, titleCornA, "x,y,z\nMode: min"], mkTri[dCornA, 6, titleCornA, "x,y,z\nMode: median"]},
     {mkTri[dCornA, 7, titleCornA, "Bernoulli\nMode: min"], mkTri[dCornA, 8, titleCornA, "Bernoulli\nMode: median"]},
     {mkTri[dCovA, 5, titleCovA, "x,y,z\nMode: min"], mkTri[dCovA, 6, titleCovA, "x,y,z\nMode: median"]},
     {mkTri[dCovA, 7, titleCovA, "Bernoulli\nMode: min"], mkTri[dCovA, 8, titleCovA, "Bernoulli\nMode: median"]}
     }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
Export[FileNameJoin[{figDir, "fig_weight2_N" <> ToString[mainOrd] <> "_" <> cornerMode <> ".png"}], checkA, ImageResolution -> figDPI];
checkA

(* ---- check B: weight 2 vs exact, evolution across w2Orders --------------- *)
evolDataCorn = (# -> load[cornerMode, 2, #, #]) & /@ w2Orders // Association;
evolDataCov = (# -> load["coverage", 2, #, #]) & /@ w2Orders // Association;
(* three separate files, same split as check C: x,y,z (cornerMode), Bernoulli
   (cornerMode), coverage (both representations together). *)
checkBxyz = Column[{Style["weight 2  vs exact  --  " <> cornerLabel <> "  x,y,z  --  convergence vs max_power = " <>
    ToString[w2Orders] <> "  (min & median)", 14, Bold],
   withLegend[Grid[{
     Table[mkTri[evolDataCorn[o], colOf[True, "min"], cornerLabel <> ", max_pow=" <> ToString[o], "x,y,z\nMode: min"], {o, w2Orders}],
     Table[mkTri[evolDataCorn[o], colOf[True, "median"], cornerLabel <> ", max_pow=" <> ToString[o], "x,y,z\nMode: median"], {o, w2Orders}]
     }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
Export[FileNameJoin[{figDir, "fig_weight2_evolution_" <> cornerMode <> "_xyz.png"}], checkBxyz, ImageResolution -> figDPI];

checkBber = Column[{Style["weight 2  vs exact  --  " <> cornerLabel <> "  Bernoulli  --  convergence vs max_power = " <>
    ToString[w2Orders] <> "  (min & median)", 14, Bold],
   withLegend[Grid[{
     Table[mkTri[evolDataCorn[o], colOf[False, "min"], cornerLabel <> ", max_pow=" <> ToString[o], "Bernoulli\nMode: min"], {o, w2Orders}],
     Table[mkTri[evolDataCorn[o], colOf[False, "median"], cornerLabel <> ", max_pow=" <> ToString[o], "Bernoulli\nMode: median"], {o, w2Orders}]
     }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
Export[FileNameJoin[{figDir, "fig_weight2_evolution_" <> cornerMode <> "_bernoulli.png"}], checkBber, ImageResolution -> figDPI];

checkBcovXyz = Column[{Style["weight 2  vs exact  --  coverage  x,y,z  --  convergence vs max_power = " <>
    ToString[w2Orders] <> "  (min & median)", 14, Bold],
   withLegend[Grid[{
     Table[mkTri[evolDataCov[o], colOf[True, "min"], "coverage, max_pow=" <> ToString[o], "x,y,z\nMode: min"], {o, w2Orders}],
     Table[mkTri[evolDataCov[o], colOf[True, "median"], "coverage, max_pow=" <> ToString[o], "x,y,z\nMode: median"], {o, w2Orders}]
     }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
Export[FileNameJoin[{figDir, "fig_weight2_evolution_" <> cornerMode <> "_coverage_xyz.png"}], checkBcovXyz, ImageResolution -> figDPI];

checkBcovBer = Column[{Style["weight 2  vs exact  --  coverage  Bernoulli  --  convergence vs max_power = " <>
    ToString[w2Orders] <> "  (min & median)", 14, Bold],
   withLegend[Grid[{
     Table[mkTri[evolDataCov[o], colOf[False, "min"], "coverage, max_pow=" <> ToString[o], "Bernoulli\nMode: min"], {o, w2Orders}],
     Table[mkTri[evolDataCov[o], colOf[False, "median"], "coverage, max_pow=" <> ToString[o], "Bernoulli\nMode: median"], {o, w2Orders}]
     }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
Export[FileNameJoin[{figDir, "fig_weight2_evolution_" <> cornerMode <> "_coverage_bernoulli.png"}], checkBcovBer, ImageResolution -> figDPI];

{checkBxyz, checkBber, checkBcovXyz, checkBcovBer}

(* ---- check C: weight 6 self-convergence, adjacent-order pairs w6Pairs ----
   Only plots pairs where BOTH cornerMode and coverage compare_*.m files
   actually exist -- run_analyze.py may not have been run for every pair yet
   (e.g. 19->20 still pending while 18->19 is done), so this degrades
   gracefully (fewer columns) instead of crashing on a missing file, and
   automatically picks up a pair once its analysis completes, no edit needed. *)
w6PairsAvail = Select[w6Pairs, haveData[cornerMode, 6, #[[1]], #[[2]]] && haveData["coverage", 6, #[[1]], #[[2]]] &];
w6PairsMissing = Complement[w6Pairs, w6PairsAvail];
If[w6PairsMissing =!= {},
   Print["[plot_compare] weight 6: skipping pair(s) with no compare_*.m yet (for ", cornerMode,
         " and/or coverage): ", w6PairsMissing]];
If[w6PairsAvail === {},
   Print["[plot_compare] weight 6: NO pairs available yet for ", cornerMode, "/coverage -- check C skipped"],
 (
  pairDataCorn = (# -> load[cornerMode, 6, #[[1]], #[[2]]]) & /@ w6PairsAvail // Association;
  pairDataCov = (# -> load["coverage", 6, #[[1]], #[[2]]]) & /@ w6PairsAvail // Association;
  pairLbl[p_] := "max_power=" <> ToString[p[[1]]] <> "\[Rule]" <> ToString[p[[2]]];
  (* NOTE: names below avoid underscores (checkCxyz not checkC_xyz) -- a bare
     name_word is ALWAYS parsed as Pattern[name, Blank[word]] by Mathematica,
     never as a plain identifier; using checkC_xyz here silently broke the
     Export (wrote out the unevaluated pattern's print-name as a tiny image
     instead of the real graphics) until caught and fixed. *)

  (* three separate files instead of one combined grid: x,y,z (cornerMode),
     Bernoulli (cornerMode), and coverage (both representations together,
     since coverage is its own separate check, not split by representation).
     Title = corner + max_pow (per column); box = representation + Mode. *)
  checkCxyz = Column[{Style["weight 6  self-convergence  --  " <> cornerLabel <> "  x,y,z  pairs " <>
      ToString[w6PairsAvail], 14, Bold],
     withLegend[Grid[{
       Table[mkTri[pairDataCorn[p], colOf[True, "min"], cornerLabel <> ", " <> pairLbl[p], "x,y,z\nMode: min"], {p, w6PairsAvail}],
       Table[mkTri[pairDataCorn[p], colOf[True, "median"], cornerLabel <> ", " <> pairLbl[p], "x,y,z\nMode: median"], {p, w6PairsAvail}]
       }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
  Export[FileNameJoin[{figDir, "fig_weight6_selfconv_" <> cornerMode <> "_xyz.png"}], checkCxyz, ImageResolution -> figDPI];

  checkCber = Column[{Style["weight 6  self-convergence  --  " <> cornerLabel <> "  Bernoulli  pairs " <>
      ToString[w6PairsAvail], 14, Bold],
     withLegend[Grid[{
       Table[mkTri[pairDataCorn[p], colOf[False, "min"], cornerLabel <> ", " <> pairLbl[p], "Bernoulli\nMode: min"], {p, w6PairsAvail}],
       Table[mkTri[pairDataCorn[p], colOf[False, "median"], cornerLabel <> ", " <> pairLbl[p], "Bernoulli\nMode: median"], {p, w6PairsAvail}]
       }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
  Export[FileNameJoin[{figDir, "fig_weight6_selfconv_" <> cornerMode <> "_bernoulli.png"}], checkCber, ImageResolution -> figDPI];

  checkCcovXyz = Column[{Style["weight 6  self-convergence  --  coverage  x,y,z  pairs " <>
      ToString[w6PairsAvail], 14, Bold],
     withLegend[Grid[{
       Table[mkTri[pairDataCov[p], colOf[True, "min"], "coverage, " <> pairLbl[p], "x,y,z\nMode: min"], {p, w6PairsAvail}],
       Table[mkTri[pairDataCov[p], colOf[True, "median"], "coverage, " <> pairLbl[p], "x,y,z\nMode: median"], {p, w6PairsAvail}]
       }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
  Export[FileNameJoin[{figDir, "fig_weight6_selfconv_" <> cornerMode <> "_coverage_xyz.png"}], checkCcovXyz, ImageResolution -> figDPI];

  checkCcovBer = Column[{Style["weight 6  self-convergence  --  coverage  Bernoulli  pairs " <>
      ToString[w6PairsAvail], 14, Bold],
     withLegend[Grid[{
       Table[mkTri[pairDataCov[p], colOf[False, "min"], "coverage, " <> pairLbl[p], "Bernoulli\nMode: min"], {p, w6PairsAvail}],
       Table[mkTri[pairDataCov[p], colOf[False, "median"], "coverage, " <> pairLbl[p], "Bernoulli\nMode: median"], {p, w6PairsAvail}]
       }, Spacings -> {1, 1}]]}, Spacings -> 1.5];
  Export[FileNameJoin[{figDir, "fig_weight6_selfconv_" <> cornerMode <> "_coverage_bernoulli.png"}], checkCcovBer, ImageResolution -> figDPI];

  {checkCxyz, checkCber, checkCcovXyz, checkCcovBer}
 )]

(* ---- check D: weight 6 improving-vs-worsening -- sign of (min_19->20 - min_18->19) ----
   Green = higher order pair has MORE agreeing digits at that point (still converging
   as N grows); red = FEWER agreeing digits (getting worse with N -- a red point is a
   real warning sign, not just "not yet converged"). MIN only (the worst-component
   digit count -- the number you'd actually trust), both representations, cornerMode
   and coverage. Requires the first two entries of w6PairsAvail (i.e. both 18->19 and
   19->20 present); degrades to a skip message otherwise, same spirit as check C. *)
If[Length[w6PairsAvail] < 2,
   Print["[plot_compare] weight 6 delta check: need >=2 available pairs, only have ", w6PairsAvail, " -- skipped"],
 (
  pLo = w6PairsAvail[[1]]; pHi = w6PairsAvail[[2]];
  cfDelta[q_] := Which[q > 0, RGBColor[0.2, 0.6, 0.25], q < 0, RGBColor[0.75, 0.15, 0.1], True, GrayLevel[0.85]];
  deltaData[d1_, d2_, col_] := Table[{d1[[i, 1]], d1[[i, 2]], d2[[i, col]] - d1[[i, col]]}, {i, Length[d1]}];
  mkTriDelta[dd_, lbl_] := Module[{hw = gridStep[dd]/2*1.02},
     Graphics[Prepend[
        Map[{cfDelta[#[[3]]], Rectangle[{#[[1]], #[[2]]} - hw, {#[[1]], #[[2]]} + hw]} &, dd],
        EdgeForm[None]],
     Frame -> True, PlotLabel -> lbl, AspectRatio -> 1, ImageSize -> 260,
     PlotRange -> {{0, 1}, {0, 1}}, FrameLabel -> {"x", "y"},
     Epilog -> {GrayLevel[0.45], Line[{{0, 0}, {1, 0}, {0, 1}, {0, 0}}]}]];
  deltaLegend = Row[{Graphics[{RGBColor[0.2, 0.6, 0.25], Rectangle[]}, ImageSize -> 14], " improving   ",
     Graphics[{RGBColor[0.75, 0.15, 0.1], Rectangle[]}, ImageSize -> 14], " worsening   ",
     Graphics[{GrayLevel[0.85], Rectangle[]}, ImageSize -> 14], " no change"}];

  ddCornXyz = deltaData[pairDataCorn[pLo], pairDataCorn[pHi], 5];
  ddCornBer = deltaData[pairDataCorn[pLo], pairDataCorn[pHi], 7];
  ddCovXyz = deltaData[pairDataCov[pLo], pairDataCov[pHi], 5];
  ddCovBer = deltaData[pairDataCov[pLo], pairDataCov[pHi], 7];

  checkDxyz = Column[{Style["weight 6  min-digit CHANGE  --  " <> cornerLabel <> "  x,y,z  --  (" <>
      pairLbl[pHi] <> ") minus (" <> pairLbl[pLo] <> ")", 14, Bold], deltaLegend,
     mkTriDelta[ddCornXyz, cornerLabel <> "  x,y,z"]}, Spacings -> 1];
  Export[FileNameJoin[{figDir, "fig_weight6_delta_" <> cornerMode <> "_xyz.png"}], checkDxyz, ImageResolution -> figDPI];

  checkDber = Column[{Style["weight 6  min-digit CHANGE  --  " <> cornerLabel <> "  Bernoulli  --  (" <>
      pairLbl[pHi] <> ") minus (" <> pairLbl[pLo] <> ")", 14, Bold], deltaLegend,
     mkTriDelta[ddCornBer, cornerLabel <> "  Bernoulli"]}, Spacings -> 1];
  Export[FileNameJoin[{figDir, "fig_weight6_delta_" <> cornerMode <> "_bernoulli.png"}], checkDber, ImageResolution -> figDPI];

  checkDcovXyz = Column[{Style["weight 6  min-digit CHANGE  --  coverage  x,y,z  --  (" <>
      pairLbl[pHi] <> ") minus (" <> pairLbl[pLo] <> ")", 14, Bold], deltaLegend,
     mkTriDelta[ddCovXyz, "coverage  x,y,z"]}, Spacings -> 1];
  Export[FileNameJoin[{figDir, "fig_weight6_delta_" <> cornerMode <> "_coverage_xyz.png"}], checkDcovXyz, ImageResolution -> figDPI];

  checkDcovBer = Column[{Style["weight 6  min-digit CHANGE  --  coverage  Bernoulli  --  (" <>
      pairLbl[pHi] <> ") minus (" <> pairLbl[pLo] <> ")", 14, Bold], deltaLegend,
     mkTriDelta[ddCovBer, "coverage  Bernoulli"]}, Spacings -> 1];
  Export[FileNameJoin[{figDir, "fig_weight6_delta_" <> cornerMode <> "_coverage_bernoulli.png"}], checkDcovBer, ImageResolution -> figDPI];

  {checkDxyz, checkDber, checkDcovXyz, checkDcovBer}
 )]


