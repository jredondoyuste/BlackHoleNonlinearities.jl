#!/usr/bin/env wolframscript

(*
  Independent Mathematica reference generator for the regularized 2nd-order
  source terms defined in mma/CoefficientsXYZ.wl (XYZ in EEO, EOE, EOO, OOE, OOO).

  Usage (run from the BlackHoleNonlinearities.jl package directory):
    wolframscript -file tools/mma_reference.wls              (* all 5 sectors *)
    wolframscript -file tools/mma_reference.wls EOE OOO       (* only these sectors *)

  Writes, per sector XYZ:
    test/reference/mma_source_XYZ.dat  -- one row per (multipole case, r, theta)
    test/reference/mma_coeffs_XYZ.dat  -- one row per multipole case (c[i,j], c[i][log], zpQ/zmQ)

  Each RegSourceXYZ expression is not edited; this script only reads
  mma/CoefficientsXYZ.wl (never modified) and substitutes numeric values.
*)

baseDir = If[$InputFileName =!= "", DirectoryName[$InputFileName, 2], Directory[]];
If[baseDir === "" || !DirectoryQ[baseDir], baseDir = Directory[]];
mmaDir = FileNameJoin[{baseDir, "mma"}];
refDir = FileNameJoin[{baseDir, "test", "reference"}];
If[!DirectoryQ[refDir], CreateDirectory[refDir, CreateIntermediateDirectories -> True]];

Print["baseDir = ", baseDir];
Print["mmaDir  = ", mmaDir];
Print["refDir  = ", refDir];

(* ---------------------------------------------------------------------- *)
(* Fixed numeric inputs (exact rationals -> high precision reference)     *)
(* ---------------------------------------------------------------------- *)

w1val = 37/100;
w2val = 52/100 - (9/100)*I;

phi1val = 42/100 + (61/100)*I;
psi1val = -37/100 + (28/100)*I;
phi2val = 85/100 - (19/100)*I;
psi2val = -53/100 - (44/100)*I;

c1m2val = 3/10 + (1/5)*I;
c1m1val = -9/20 + (13/20)*I;
c2m2val = 11/20 - (7/20)*I;
c2m1val = -1/4 - (1/2)*I;
c3m2val = 2/5 + (3/20)*I;
c4m2val = -3/5 + (1/10)*I;

thetaVals = {7/10, 13/10};
rVals = {5, 10};

(* per-sector override of rVals, used only if EEO turns out too slow *)
sectorRValsOverride = <||>;

(* ---------------------------------------------------------------------- *)
(* Formatting helpers -- 17 significant digit scientific notation,        *)
(* readable by Julia's parse(Float64, ...) / readdlm.                     *)
(* ---------------------------------------------------------------------- *)

FormatReal[x_?NumericQ] := Module[{n, sign, absn, exp, mant, intPart, digitsStr},
  n = N[x, 25];
  If[n == 0, Return["0.0000000000000000e+00"]];
  sign = If[n < 0, "-", ""];
  absn = Abs[n];
  exp = Floor[Log10[absn]];
  mant = absn/10^exp;
  If[mant >= 10, mant = mant/10; exp = exp + 1];
  If[mant < 1, mant = mant*10; exp = exp - 1];
  intPart = Round[mant*10^16];
  If[intPart >= 10^17, intPart = Round[intPart/10]; exp = exp + 1];
  digitsStr = IntegerString[intPart, 10, 17];
  StringJoin[sign, StringTake[digitsStr, 1], ".", StringDrop[digitsStr, 1],
    "e", If[exp >= 0, "+", "-"], IntegerString[Abs[exp], 10, 2]]
]

RowJoin[fields_List] := StringJoin[Riffle[fields, " "]]

SafeThreeJSymbol[a_, b_, cc_, d_, e_, f_] := Module[{v},
  v = Quiet[ThreeJSymbol[{a, b}, {cc, d}, {e, f}]];
  If[NumericQ[v], v, 0]
]

jsRule = JS[{a_, b_}, {cc_, d_}, {e_, f_}] :> SafeThreeJSymbol[a, b, cc, d, e, f];

(* flag anything whose tracked precision fell well short of the requested
   30 digits -- catches genuine cancellation problems that Off[N::meprec]
   would otherwise hide *)
LowPrecQ[v_] := !NumericQ[v] || (Precision[v] =!= Infinity && Precision[v] < 12)

(* Evaluate substitutedExpr to 30 digits; if the tracked precision is too
   low (deep radical cancellation -- verified on spot checks to sometimes
   be an exact 0 that N[] alone can't resolve cleanly), retry once with
   Simplify first (time-boxed: the c[i,j]/zpQ/zmQ formulas are small so
   this is fast, but a fallback exists in case it's ever used on something
   larger). Returns {value, retried, stillBad}. *)
RobustN[substitutedExpr_] := Module[{val, val2},
  val = N[substitutedExpr, 30];
  If[LowPrecQ[val],
    val2 = TimeConstrained[Quiet[N[Simplify[substitutedExpr], 30]], 60, $Failed];
    If[NumericQ[val2] && !LowPrecQ[val2],
      {val2, True, False},
      {val, True, True}
    ],
    {val, False, False}
  ]
]

(* the c[i,j]/zpQ/zmQ formulas contain many nested Sqrt[] terms with large
   cancellations; the default $MaxExtraPrecision (50) is not enough for
   N[...,30] to converge cleanly on some of them, raising N::meprec.
   Raise the cap so the requested 30 digits are genuinely reliable. *)
$MaxExtraPrecision = 2000;
(* Even with $MaxExtraPrecision raised, N::meprec still fires on some
   deeply-nested-radical c[i,j]/zpQ/zmQ formulas; spot checks (see
   tools/ dev notes) show the returned value still carries ~30 digits of
   Precision[] regardless, so the message is cosmetic here -- silence it. *)
Off[N::meprec];

(* kept EXACT (no premature N[]) so the final N[...,30] on the fully
   substituted expression doesn't hit mixed exact/approximate precision
   loss when combined with the many nested Sqrt[] factors in the c[i,j] *)
KKvalExact[L1_, M1_, L2_, M2_, L_] :=
  SafeThreeJSymbol[L1, M1, L2, M2, L, -(M1 + M2)]*
     Sqrt[(2*L + 1)*(2*L1 + 1)*(2*L2 + 1)/(4*Pi)];

(* ---------------------------------------------------------------------- *)
(* Sector configuration                                                    *)
(* ---------------------------------------------------------------------- *)

oddSumCases  = {{2, 2, 2, -2, 3}, {2, 1, 3, 0, 4}, {3, 2, 2, -1, 2}};
evenSumCases = {{2, 2, 2, -2, 2}, {2, 1, 3, 0, 3}, {3, 2, 2, -1, 3}, {2, 2, 2, 2, 4}};

sectorConfig = <|
  "EEO" -> <|"file" -> "CoefficientsEEO.wl", "reg" -> "RegSourceEEO",
             "field1" -> "\[Psi]1l1e", "field2" -> "\[Psi]1l2e",
             "cases" -> oddSumCases, "zsym" -> "zmQ"|>,
  "EOE" -> <|"file" -> "CoefficientsEOE.wl", "reg" -> "RegSourceEOE",
             "field1" -> "\[Psi]1l1e", "field2" -> "\[Psi]1l2o",
             "cases" -> oddSumCases, "zsym" -> "zpQ"|>,
  "EOO" -> <|"file" -> "CoefficientsEOO.wl", "reg" -> "RegSourceEOO",
             "field1" -> "\[Psi]1l1e", "field2" -> "\[Psi]1l2o",
             "cases" -> evenSumCases, "zsym" -> "zmQ"|>,
  "OOE" -> <|"file" -> "CoefficientsOOE.wl", "reg" -> "RegSourceOOE",
             "field1" -> "\[Psi]1l1o", "field2" -> "\[Psi]1l2o",
             "cases" -> evenSumCases, "zsym" -> "zpQ"|>,
  "OOO" -> <|"file" -> "CoefficientsOOO.wl", "reg" -> "RegSourceOOO",
             "field1" -> "\[Psi]1l1o", "field2" -> "\[Psi]1l2o",
             "cases" -> oddSumCases, "zsym" -> "zmQ"|>
|>;

cIndexName[j_Integer] := If[j < 0, "m" <> ToString[-j], ToString[j]];

(* ---------------------------------------------------------------------- *)
(* Per-sector processing. All file-defined symbols are localized with     *)
(* Block so successive Get[] calls never leak definitions between sectors *)
(* (equivalent in effect to the ClearAll["Global`*"] the task describes,  *)
(* but without also wiping this driver script's own variables).          *)
(* ---------------------------------------------------------------------- *)

ProcessSector[label_String] := Module[
  {cfg, sourcePath, coeffPath, sourceStream, coeffStream, t0, t1, rValsHere},

  cfg = sectorConfig[label];
  rValsHere = Lookup[sectorRValsOverride, label, rVals];

  sourcePath = FileNameJoin[{refDir, "mma_source_" <> label <> ".dat"}];
  coeffPath  = FileNameJoin[{refDir, "mma_coeffs_" <> label <> ".dat"}];

  sourceStream = OpenWrite[sourcePath];
  WriteLine[sourceStream,
    "# l1 m1 l2 m2 l re_w1 im_w1 re_w2 im_w2 theta r re_phi1 im_phi1 re_psi1 im_psi1 re_phi2 im_phi2 re_psi2 im_psi2 re_c1m2 im_c1m2 re_c1m1 im_c1m1 re_c2m2 im_c2m2 re_c2m1 im_c2m1 re_c3m2 im_c3m2 re_c4m2 im_c4m2 re_S im_S"];

  coeffStream = OpenWrite[coeffPath];

  t0 = AbsoluteTime[];

  Block[{RegSourceEEO, RegSourceEOE, RegSourceEOO, RegSourceOOE, RegSourceOOO,
         c, JS, KK, zpQ, zmQ, \[Lambda]2l, parity, M, l, l1, l2, r,
         \[Omega]1, \[Omega]2, \[Theta], \[Chi]1, \[Chi]2,
         \[Psi]1l1e, \[Psi]1l1o, \[Psi]1l2e, \[Psi]1l2o},

    Module[
      {regSym, field1Sym, field2Sym, zSym, expr, cPairs, cLogIdx, zTriples,
       coeffCols, coeffHeaderFields, warnCount = 0},

      Get[FileNameJoin[{mmaDir, cfg["file"]}]];

      regSym    = Symbol[cfg["reg"]];
      field1Sym = Symbol[cfg["field1"]];
      field2Sym = Symbol[cfg["field2"]];
      zSym      = Symbol[cfg["zsym"]];

      expr = regSym; (* forces full expansion of the c[...] definitions *)

      (* NB: DownValues/SubValues store LHS's wrapped in a literal HoldPattern
         head. Plain HoldPattern[...] used *inside* a Cases pattern is
         transparent (it matches what its argument would match), so it does
         NOT match that literal wrapper -- Verbatim[HoldPattern] is required.
         zSym is stored in a variable, so DownValues needs Evaluate[] to see
         through to the actual symbol (DownValues has HoldFirst). *)
      cPairs   = Cases[DownValues[c], (Verbatim[HoldPattern][c[i_Integer, j_Integer]] :> _) :> {i, j}];
      cLogIdx  = Cases[SubValues[c], (Verbatim[HoldPattern][_[i_Integer][log]] :> _) :> i];
      zTriples = Cases[DownValues[Evaluate[zSym]], (Verbatim[HoldPattern][zSym[i_Integer, j_Integer, k_Integer]] :> _) :> {i, j, k}];

      Print[label, ": ", Length[cPairs], " defined c[i,j], ", Length[cLogIdx],
        " defined c[i][log], ", Length[zTriples], " defined ", cfg["zsym"], "[i,j,k]"];

      coeffCols = Join[
        ("c_" <> ToString[#[[1]]] <> "_" <> cIndexName[#[[2]]]) & /@ cPairs,
        ("clog_" <> ToString[#]) & /@ cLogIdx,
        ("zq_" <> ToString[#[[1]]] <> "_" <> ToString[#[[2]]] <> "_" <> cIndexName[#[[3]]]) & /@ zTriples
      ];
      coeffHeaderFields = Join[{"l1", "m1", "l2", "m2", "l", "re_w1", "im_w1", "re_w2", "im_w2"},
        Flatten[{"re_" <> #, "im_" <> #} & /@ coeffCols]];
      WriteLine[coeffStream, "# " <> RowJoin[coeffHeaderFields]];

      Do[
        Module[{L1, M1, L2, M2, L, kk, multRules, cFreeRules, coeffTokens, val, coeffExpr, retried, stillBad},
          {L1, M1, L2, M2, L} = testcase;
          multRules = {l1 -> L1, l2 -> L2, l -> L};
          kk = KKvalExact[L1, M1, L2, M2, L];
          cFreeRules = {c[1, -2] -> c1m2val, c[1, -1] -> c1m1val, c[2, -2] -> c2m2val,
                        c[2, -1] -> c2m1val, c[3, -2] -> c3m2val, c[4, -2] -> c4m2val};

          (* ---- coefficients row ---- *)
          coeffTokens = {ToString[L1], ToString[M1], ToString[L2], ToString[M2], ToString[L],
            FormatReal[Re[w1val]], FormatReal[Im[w1val]], FormatReal[Re[w2val]], FormatReal[Im[w2val]]};

          Do[
            coeffExpr = c[pair[[1]], pair[[2]]];
            {val, retried, stillBad} = RobustN[coeffExpr /. cFreeRules /. multRules /. KK -> kk /.
                     {\[Omega]1 -> w1val, \[Omega]2 -> w2val} /. jsRule];
            If[stillBad, warnCount++;
              Print["  WARNING [", label, "] low-precision/non-numeric c[", pair[[1]], ",", pair[[2]],
                    "] at case ", testcase, ": Precision=", If[NumericQ[val], Precision[val], "N/A"]]];
            coeffTokens = Join[coeffTokens, {FormatReal[Re[val]], FormatReal[Im[val]]}];
            , {pair, cPairs}];

          Do[
            coeffExpr = c[idx][log];
            {val, retried, stillBad} = RobustN[coeffExpr /. cFreeRules /. multRules /. KK -> kk /.
                     {\[Omega]1 -> w1val, \[Omega]2 -> w2val} /. jsRule];
            If[stillBad, warnCount++;
              Print["  WARNING [", label, "] low-precision/non-numeric c[", idx, "][log] at case ", testcase,
                    ": Precision=", If[NumericQ[val], Precision[val], "N/A"]]];
            coeffTokens = Join[coeffTokens, {FormatReal[Re[val]], FormatReal[Im[val]]}];
            , {idx, cLogIdx}];

          Do[
            coeffExpr = zSym[trip[[1]], trip[[2]], trip[[3]]];
            {val, retried, stillBad} = RobustN[coeffExpr /. cFreeRules /. multRules /. KK -> kk /.
                     {\[Omega]1 -> w1val, \[Omega]2 -> w2val, \[Chi]1 -> 1, \[Chi]2 -> -1} /. jsRule];
            If[stillBad, warnCount++;
              Print["  WARNING [", label, "] low-precision/non-numeric ", cfg["zsym"], "[", trip[[1]], ",", trip[[2]], ",", trip[[3]],
                    "] at case ", testcase, ": Precision=", If[NumericQ[val], Precision[val], "N/A"]]];
            coeffTokens = Join[coeffTokens, {FormatReal[Re[val]], FormatReal[Im[val]]}];
            , {trip, zTriples}];

          WriteLine[coeffStream, RowJoin[coeffTokens]];

          (* ---- source rows, over r and theta ---- *)
          Do[
            Do[
              Module[{fieldRules, e1, sVal, tokens, sRetried, sStillBad},
                fieldRules = {
                  Derivative[1][field1Sym][r] -> psi1val, field1Sym[r] -> phi1val,
                  Derivative[1][field2Sym][r] -> psi2val, field2Sym[r] -> phi2val
                };
                e1 = expr /. fieldRules;
                e1 = e1 /. cFreeRules;
                e1 = e1 /. multRules;
                e1 = e1 /. KK -> kk;
                e1 = e1 /. {\[Omega]1 -> w1val, \[Omega]2 -> w2val};
                e1 = e1 /. \[Theta] -> thV;
                e1 = e1 /. r -> rV;
                e1 = e1 /. jsRule;
                {sVal, sRetried, sStillBad} = RobustN[e1];
                If[sStillBad, warnCount++;
                  Print["  WARNING [", label, "] low-precision/non-numeric source at case ", testcase,
                        " r=", rV, " theta=", thV, ": Precision=", If[NumericQ[sVal], Precision[sVal], "N/A"]]];

                tokens = {ToString[L1], ToString[M1], ToString[L2], ToString[M2], ToString[L],
                  FormatReal[Re[w1val]], FormatReal[Im[w1val]], FormatReal[Re[w2val]], FormatReal[Im[w2val]],
                  FormatReal[thV], FormatReal[rV],
                  FormatReal[Re[phi1val]], FormatReal[Im[phi1val]],
                  FormatReal[Re[psi1val]], FormatReal[Im[psi1val]],
                  FormatReal[Re[phi2val]], FormatReal[Im[phi2val]],
                  FormatReal[Re[psi2val]], FormatReal[Im[psi2val]],
                  FormatReal[Re[c1m2val]], FormatReal[Im[c1m2val]],
                  FormatReal[Re[c1m1val]], FormatReal[Im[c1m1val]],
                  FormatReal[Re[c2m2val]], FormatReal[Im[c2m2val]],
                  FormatReal[Re[c2m1val]], FormatReal[Im[c2m1val]],
                  FormatReal[Re[c3m2val]], FormatReal[Im[c3m2val]],
                  FormatReal[Re[c4m2val]], FormatReal[Im[c4m2val]],
                  FormatReal[Re[sVal]], FormatReal[Im[sVal]]};
                WriteLine[sourceStream, RowJoin[tokens]];
              ]
              , {thV, thetaVals}]
            , {rV, rValsHere}]
        ]
        , {testcase, cfg["cases"]}];

      If[warnCount > 0, Print["  ", label, ": ", warnCount, " WARNING(S) -- check output for non-numeric / unevaluated symbols"]];
    ]
  ]; (* end Block *)

  Close[sourceStream]; Close[coeffStream];
  t1 = AbsoluteTime[];
  Print[label, " done in ", N[t1 - t0, 6], " s  -> ", sourcePath, "  ", coeffPath];
  t1 - t0
]

(* ---------------------------------------------------------------------- *)
(* Driver                                                                  *)
(* ---------------------------------------------------------------------- *)

argSectors = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], {}];
allSectors = {"EEO", "EOE", "EOO", "OOE", "OOO"};
sectorsToRun = If[argSectors === {}, allSectors, Intersection[allSectors, ToUpperCase /@ argSectors]];

Print["Sectors to run: ", sectorsToRun];

timings = <||>;
Do[
  timings[label] = ProcessSector[label];
  , {label, sectorsToRun}];

Print["---- timings (s) ----"];
Print[timings];
