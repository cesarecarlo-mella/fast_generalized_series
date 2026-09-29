(* ============================================================================
   00_transport_boundaries.wl  --  boundary constants at CC1/CC2/CC3, family H
   ----------------------------------------------------------------------------
   CC1, CC2, CC3 are three NEW kinematic regions, all reached by analytically
   continuing OUT from the EXISTING corner-1 boundary data -- they are NOT the
   analytic continuation of the existing blow-up corners C2/C3 used by
   ../../../all_families/ (careful: same "2"/"3" labels, unrelated regions).
   Output uses "CC" (double-C) throughout specifically to avoid that confusion
   with all_families/<Family>/pipeline/dist/bdy_c<c>_w<w>.m.

   CC1 -- "easy" direct continuation from corner 1 itself: take the EXISTING
          phys_c1_w<w>_ord20.m (already-computed corner-1 physical series),
          abstract Log[x]/Log[y] as placeholders, re-expand them through the
          analytic-continuation substitution LogX -> Log[w]-Log[1-w-u]+I*Pi
          (same idea for LogY), Taylor-expand the WHOLE physical series in the
          new (u,w) around 0, drop whatever logs survive, read off the
          constant. No atilde/DEQ/GIntegrate needed at all for this one.

   CC2, CC3 -- reached via the "high energy" (z->infinity) reparametrization
          x -> -w/(1-u-w), y -> -u/(1-u-w): derive the effective 1D DEQ along
          w (u=0 fixed) and along u (w=0 fixed) from the ORIGINAL 20-letter
          alphabet, dot into atilde, GIntegrate weight-by-weight (seeded at
          each weight by CC1's OWN value -- w=0,u=0 in this parametrization
          IS corner 1), then land at the path's far endpoint via the SAME
          change-of-fibration-basis + regularization used throughout this
          project: reexpress as 1-t (safety), t->0, throw away the resulting
          divergent all-zero-string G[0,...,0,t] (finite otherwise -- verified
          weight-by-weight, all 126 fibrtto1.m rules through weight 6, no
          leftover G[...] or t-dependence survives this regularized limit).

   Both routines (pathW/pathU construction, and fibrtto1.m's regularization
   completeness) were independently cross-checked before writing this script:
     - fibrtto1.m: all 126 rules (weight 1..6, complete {0,1}-alphabet)
       regularize-and-limit to a pure constant, verified via GetGs.
     - pathW/pathU: 18/20 letters match a naive independent differentiation
       exactly; the 2 sqrt-sector exceptions (l7,l8) were verified via an
       actual numerical limit (not just the symbolic Series shortcut).

   OUTPUT (into ./dist/, SAME convention as all_families/.../pipeline/dist/):
       bdy_CC1_w<w>.m , bdy_CC2_w<w>.m , bdy_CC3_w<w>.m   w = 0..6
       each a length-Nf list of pure Pi/Zeta constants (like bdy_c<c>_w<w>.m)

   USAGE :  wolframscript -file 00_transport_boundaries.wl  [TOWEIGHT]
            default TOWEIGHT=6; pass e.g. 1 for a cheap partial/sanity run
            (positional, NOT a --flag -- wolframscript -file silently drops
            any --flag argument instead of forwarding it to the script)
   ============================================================================ *)

Needs["PolyLogTools`"];

scriptPath = Which[
   $ScriptCommandLine =!= {}, First[$ScriptCommandLine],
   $InputFileName =!= "", $InputFileName,
   True, None];
If[scriptPath === None,
   Print["ERROR: cannot determine this script's own file path. Run via ",
         "`wolframscript -file 00_transport_boundaries.wl`."];
   Exit[1]];
famDir = DirectoryName[AbsoluteFileName[scriptPath]];             (* .../transport.boundaries/H/ *)
projectRoot = ParentDirectory[ParentDirectory[ParentDirectory[famDir]]]; (* .../series.all *)
SetDirectory[projectRoot];
dist = FileNameJoin[{famDir, "dist"}] <> "/";
Quiet@CreateDirectory[dist];

FAMILY = "H";

argv = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];
TOWEIGHT = If[Length[argv] >= 1, ToExpression[argv[[1]]], 6];

Print["=== 00_transport_boundaries : family = ", FAMILY, "  weights 1..", TOWEIGHT, " ==="];
t0all = AbsoluteTime[];

(* ---- inputs ---- *)
atilde = Import["../all.ints/a_tildes/" <> FAMILY <> "_atilde.m"];
Nf = Length[atilde];
lettersdict = Import["../all.ints/a_tildes/letters.dict.m"];
redtto1 = Import["change_of_fibrationBasis/fibrtto1.m"];
bc = Import["../all.ints/boundaries_all/" <> FAMILY <> "_bc.m"];
sol0 = Quiet[Coefficient[#, ep, 0] & /@ bc[[1 ;; Nf, 2]]];
Print["[00] Nf=", Nf, "  sol0 nonzero: ", Count[sol0, x_ /; x =!= 0], "/", Nf];

zeroRules = Table[G @@ Join[ConstantArray[0, n], {t}] -> 0, {n, 1, 6}];

(* after GIntegrate + /. redtto1, every surviving G[...] should be a pure
   function of t (redtto1's whole job is turning G[...,1-t] into G[...,t]) --
   check this explicitly rather than assuming it, since it would silently
   miss any letter/word redtto1's 126-rule table doesn't cover. *)
checkAllT[expr_, label_String] := Module[{gs, bad},
   gs = GetGs[expr];
   bad = Select[gs, Function[g, Last[List @@ g] =!= t]];
   If[bad =!= {},
      Print["  !! WARNING (", label, "): ", Length[bad],
            " G[...] NOT purely a function of t (still involve 1-t/w/u/other): ",
            ToString[Take[bad, UpTo[5]], InputForm]],
      Print["  [check] ", label, ": all ", Length[gs], " surviving G[...] confirmed pure functions of t"]];
   ];

Put[sol0, dist <> "bdy_CC1_w0.m"];
Put[sol0, dist <> "bdy_CC2_w0.m"];
Put[sol0, dist <> "bdy_CC3_w0.m"];

(* ============================================================================
   CC1 : direct continuation from corner 1 (no atilde/DEQ needed)
   ============================================================================ *)
Print["\n--- CC1 (direct continuation from corner 1) ---"];
constC1 = <|0 -> sol0|>;
Do[
 Module[{physw, tmp, t0 = AbsoluteTime[]},
  physw = Import["cluster_results/blow_up_form_2_order_20/out/phys_c1_w" <>
                  ToString[wt] <> "_ord20.m"] /.
          {x^a_ /; a > 2 -> 0, y^a_ /; a > 2 -> 0};
  tmp = Series[Series[
      physw /. Log[x] -> LogX /. Log[y] -> LogY /. {x -> 0, y -> 0} /.
        {LogX -> Log[w] - Log[1 - w - u] + I Pi,
         LogY -> Log[u] - Log[1 - w - u] + I Pi},
      {w, 0, 2}], {u, 0, 2}] // Normal // Expand;
  constC1[wt] = (tmp /. _Log -> 0 /. {u -> 0, w -> 0}) // Expand;
  Put[constC1[wt], dist <> "bdy_CC1_w" <> ToString[wt] <> ".m"];
  Print["  weight ", wt, "  nonzero: ", Count[constC1[wt], x_ /; x =!= 0], "/", Nf,
        "  (", Round[AbsoluteTime[] - t0, .1], " s)"];
 ],
 {wt, 1, TOWEIGHT}];
Print["[00] CC1 written  w0..w", TOWEIGHT];

(* ============================================================================
   pathW / pathU : effective 1D transport connections toward CC2 / CC3
   ============================================================================ *)
Print["\n--- building pathW/pathU (effective transport connections) ---"];
tP0 = AbsoluteTime[];
subs = {x -> -w/(1 - u - w), y -> -u/(1 - u - w)};

tmpW = PowerExpand[Series[lettersdict[[;; , 2]] /. subs, {u, 0, 0},
    Assumptions -> (Im[u w/(-1 + w)] +
        Re[u w/((-1 + w)^2 Sqrt[-u w/(-1 + w)^3])] > 0)] // Normal];
pathW = Thread[Rule[lettersdict[[;; , 1]], D[tmpW, w]]];

tmpU = PowerExpand[Series[lettersdict[[;; , 2]] /. subs, {w, 0, 0},
    Assumptions -> {(Im[w/((-1 + u) u)] -
         Re[w/((-1 + u)^2 Sqrt[-u w/(-1 + u)^3])] > 0),
       (Im[u w/(-1 + u)] +
         Re[u w/((-1 + u)^2 Sqrt[-u w/(-1 + u)^3])] > 0)}] // Normal];
pathU = Thread[Rule[lettersdict[[;; , 1]], D[tmpU, u]]];

effeciteDEQw = atilde /. pathW;
effeciteDEQu = atilde /. pathU;
Print["  effeciteDEQw, effeciteDEQu built (", Nf, " x ", Nf, ")  (",
      Round[AbsoluteTime[] - tP0, .1], " s)"];

(* ============================================================================
   CC2 : transport along w -> 1 (u=0 fixed), seeded by CC1 at each weight
   ============================================================================ *)
Print["\n--- CC2 (transport along w, seeded by CC1) ---"];
constC2 = <|0 -> sol0|>;
int = sol0;
Do[
 Module[{toint, int1, afterFibration, t0 = AbsoluteTime[]},
  toint = (effeciteDEQw . int) // Expand;
  int1 = (GIntegrate[#, w] &) /@ toint;
  int = int1 + constC1[wt];
  afterFibration = (int /. w -> 1 - t) /. redtto1;
  checkAllT[afterFibration, "CC2 weight " <> ToString[wt]];
  constC2[wt] = (afterFibration /. zeroRules /. t -> 0) // Expand;
  Put[constC2[wt], dist <> "bdy_CC2_w" <> ToString[wt] <> ".m"];
  Print["  weight ", wt, "  nonzero: ", Count[constC2[wt], x_ /; x =!= 0], "/", Nf,
        "  (", Round[AbsoluteTime[] - t0, .1], " s)"];
 ],
 {wt, 1, TOWEIGHT}];
Print["[00] CC2 written  w0..w", TOWEIGHT];

(* ============================================================================
   CC3 : transport along u -> 1 (w=0 fixed), seeded by CC1 at each weight
   ============================================================================ *)
Print["\n--- CC3 (transport along u, seeded by CC1) ---"];
constC3 = <|0 -> sol0|>;
int = sol0;
Do[
 Module[{toint, int1, afterFibration, t0 = AbsoluteTime[]},
  toint = (effeciteDEQu . int) // Expand;
  int1 = (GIntegrate[#, u] &) /@ toint;
  int = int1 + constC1[wt];
  afterFibration = (int /. u -> 1 - t) /. redtto1;
  checkAllT[afterFibration, "CC3 weight " <> ToString[wt]];
  constC3[wt] = (afterFibration /. zeroRules /. t -> 0) // Expand;
  Put[constC3[wt], dist <> "bdy_CC3_w" <> ToString[wt] <> ".m"];
  Print["  weight ", wt, "  nonzero: ", Count[constC3[wt], x_ /; x =!= 0], "/", Nf,
        "  (", Round[AbsoluteTime[] - t0, .1], " s)"];
 ],
 {wt, 1, TOWEIGHT}];
Print["[00] CC3 written  w0..w", TOWEIGHT];

Print["\nDONE  (", Round[AbsoluteTime[] - t0all, .1], " s total)  -> ", dist];
