(* ============================================================================
   05_bernoulli.wl  --  STEP 5 : Bernoulli (conformal) resummation of the
   physical-variable series, per component slice.
   ----------------------------------------------------------------------------
   Runs in the CURRENT WORKING DIRECTORY.  N processes (slice p: Mod[i-1,N]==p).
   Independent of steps 1-4: reads the already-converted phys_c<c>_w<w>_ord<ORDp>.m
   (produced by 03_to_physical_vars.wl / run_to_physical.py, a series in the
   physical variables) and re-expands each component in the conformal
   variables

       x -> 1 - Exp[-sx]     (<=>  sx = -Log[1-x])
       y -> 1 - Exp[-sy]     (<=>  sy = -Log[1-y])
       z -> 1 - Exp[-sz]     (<=>  sz = -Log[1-z])

   Only the corner's OWN two physical variables are transformed (same table
   as 03_to_physical_vars.wl's SMALL/BIG):
       c1: x, y      c2: z, y      c3: x, z
   `y` is treated identically regardless of which corner it appears in --
   always y -> sy, nothing corner-specific.

   Truncated to order ORDb in EACH new s-variable present (blow_up_form(_2)
   analog of the old pipeline's 06_emit_bernoulli.wl "ORD in each of
   sIN,sOU"). Explicit Log[x],Log[y],Log[z] are protected (LogX,LogY,LogZ
   placeholders) BEFORE the substitution so the branch-cut series
   Series[Log[1-Exp[-s]],...] is never built -- this is what keeps it cheap.
   LogX/LogY/LogZ are left IN THE OUTPUT (not restored) -- a consumer
   substitutes LogX -> Log[x], LogY -> Log[y], LogZ -> Log[z] at the numeric
   point, same convention the old pipeline's Bernoulli step used.

       reads :  phys_c<c>_w<w>_ord<ORDp>.m
       writes:  ber_c<c>_w<w>_ordp<ORDp>_ordb<ORDb>_p<p>_<N>.m   { {i, r_i}, ... }
                prog05_c<c>_w<w>_p<p>.txt

   USAGE :  <kernel> 05_bernoulli.wl  <c> <w> <p> <N> <ORDp> <ORDb>
   ============================================================================ *)

Quiet[Off[End::noctx]];
cliArgs = Module[{s = $ScriptCommandLine, c = $CommandLine, pos},
   If[Length[s] > 1, Rest[s],
      pos = FirstPosition[c, f_ /; StringQ[f] && StringMatchQ[f, ___ ~~ "." ~~ ("wl"|"m")], {0}][[1]];
      If[IntegerQ[pos] && pos >= 1 && pos < Length[c], Drop[c, pos], {}]]];

If[Length[cliArgs] >= 6,
   {corner, wTarget, part, nParts, ORDp, ORDb} = ToExpression /@ Take[cliArgs, 6],
   {corner, wTarget, part, nParts, ORDp, ORDb} = {1, 1, 0, 1, 15, 15}
];

mapOf = <|1 -> {x, y}, 2 -> {z, y}, 3 -> {x, z}|>;
If[! KeyExistsQ[mapOf, corner],
   Print["[05] ERROR: corner must be 1, 2 or 3, got ", corner]; Exit[1]];
vars = mapOf[corner];
sOf = <|x -> sx, y -> sy, z -> sz|>;

cw   = "c" <> ToString[corner] <> "_w" <> ToString[wTarget];
cwpb = cw <> "_ordp" <> ToString[ORDp] <> "_ordb" <> ToString[ORDb] <>
       "_p" <> ToString[part] <> "_" <> ToString[nParts];
progFile = "prog05_c" <> ToString[corner] <> "_w" <> ToString[wTarget] <>
           "_p" <> ToString[part] <> ".txt";
Print["[05 ", cwpb, "]  vars ", vars, " -> ", sOf /@ vars, "  ORDb ", ORDb, "]  dir ", Directory[]];

inFile = "phys_" <> cw <> "_ord" <> ToString[ORDp] <> ".m";
If[! FileExistsQ[inFile], Print["[05 ", cwpb, "]  ERROR: missing ", inFile]; Exit[1]];
phys = Get[inFile];
ncomp = Length[phys];
mine = Select[Range[ncomp], Mod[# - 1, nParts] == part &];
Print["[05 ", cwpb, "]  ", Length[mine], "/", ncomp, " components in this slice"];

berify[e_] := Module[{b},
   If[PossibleZeroQ[e], Return[0]];
   b = e /. {Log[x] -> LogX, Log[y] -> LogY, Log[z] -> LogZ};
   Do[
      b = Normal@Series[b /. v -> 1 - Exp[-sOf[v]], {sOf[v], 0, ORDb}],
      {v, vars}];
   Collect[b, {LogX, LogY, LogZ, sx, sy, sz}, Together]
];

t0 = AbsoluteTime[]; out = {}; nd = 0; nn = Length[mine];
Do[
   tc = AbsoluteTime[];
   r = berify[phys[[i]]];
   dt = AbsoluteTime[] - tc;
   AppendTo[out, {i, r}];
   nd++;
   If[dt > 30, Print["[05 ", cwpb, "]  slow comp ", i, " : ", Round[dt, 0.1], " s"]];
   If[Mod[nd, 5] == 0 || nd == nn,
      Put[ToString[nd] <> " " <> ToString[nn], progFile];
      Print["[05 ", cwpb, "]  ", nd, "/", nn, "  (comp ", i, ")  ",
            Round[AbsoluteTime[] - t0, 0.1], " s"];
   ],
   {i, mine}
];
Put[out, "ber_" <> cwpb <> ".m"];
Print["[05 ", cwpb, "] DONE  ", Round[AbsoluteTime[] - t0, .1], " s"];
