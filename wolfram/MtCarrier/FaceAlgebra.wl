(* ::Package:: *)
(* Discrete Euler character C_d, S_d on Z_{N_ring(d)}: |C_d| + |S_d| == Q_sq(d) *)

IsZeroNatOp[n_Integer] := Quotient[1, Abs[n] + 1];
IsPosNatOp[n_Integer]  := 1 - IsZeroNatOp[Max[0, n]];
EqIntIndicatorOp[a_Integer, b_Integer] := IsZeroNatOp[Abs[a - b]];
BlendByIndicatorOp[sel_Integer, valTrue_, valFalse_] := sel * valTrue + (1 - sel) * valFalse;

Wrap512Op[d_Integer, phi_Integer] := Module[
  {u = Mod[phi, NRingOp[d]], q},
  q = Quotient[u + (PhiPauliOp[d] - MuVacOp[d]), NRingOp[d]];
  u - NRingOp[d] * q
];

ProjCOp[d_Integer, phi_Integer] := QSquareOp[d] - Abs[Wrap512Op[d, phi]];
ProjSOp[d_Integer, phi_Integer] := ProjCOp[d, phi - QSquareOp[d]];

(* Z_{N_ring(d)} =~ Z_4^Gray |>< Z_{Q_sq(d)} *)
GrayButterflyPairOp[d_Integer, phi_Integer] := Module[
  {u, q1, q0, m7, comp7, cMag, sMag, cNeg, sNeg, two = C0OverCSqOp[d], mu = MuVacOp[d]},
  u     = Mod[phi, NRingOp[d]];
  q1    = Quotient[u, PhiPauliOp[d]];
  q0    = Mod[Quotient[u, QSquareOp[d]], two];
  m7    = Mod[u, QSquareOp[d]];
  comp7 = QSquareOp[d] - m7;
  cMag  = (mu - q0) * comp7 + q0 * m7;
  sMag  = (mu - q0) * m7 + q0 * comp7;
  cNeg  = Mod[q1 + q0, two];
  sNeg  = q1;
  {(mu - two * cNeg) * cMag, (mu - two * sNeg) * sMag}
];

DifferenceWronskianOp[d_Integer, phi_Integer] :=
  ProjCOp[d, phi] * ProjSOp[d, phi + MuVacOp[d]] -
  ProjSOp[d, phi] * ProjCOp[d, phi + MuVacOp[d]];

PythagoreanLiftOp[d_Integer, phi_Integer] := Module[
  {c = ProjCOp[d, phi], s = ProjSOp[d, phi], q2 = QSquareOp[d]^C0OverCSqOp[d],
   s2, cNonZero, cPos, sgnC, two = C0OverCSqOp[d], mu = MuVacOp[d]},
  s2 = s^two;
  cNonZero = mu - IsZeroNatOp[Abs[c]];
  cPos = IsPosNatOp[c];
  sgnC = cNonZero * (two * cPos - mu);
  {sgnC * (q2 - s2), two * QSquareOp[d] * s, q2 + s2}
];

FresnelCosSumOp[d_Integer] :=
  Sum[ProjCOp[d, phi^C0OverCSqOp[d]], {phi, 0, NRingOp[d] - MuVacOp[d]}];

FresnelSinSumOp[d_Integer] :=
  Sum[ProjSOp[d, phi^C0OverCSqOp[d]], {phi, 0, NRingOp[d] - MuVacOp[d]}];

VerifyFaceAlgebraModuleOp[d_Integer] := Module[{ringStates, qSq = QSquareOp[d], pauli = PhiPauliOp[d]},
  ringStates = Range[0, NRingOp[d] - MuVacOp[d]];
  AllTrue[ringStates, (GrayButterflyPairOp[d, #] == {ProjCOp[d, #], ProjSOp[d, #]}) &] &&
  AllTrue[ringStates, (Abs[ProjCOp[d, #]] + Abs[ProjSOp[d, #]] == qSq) &] &&
  AllTrue[ringStates, (DifferenceWronskianOp[d, #] == qSq) &] &&
  AllTrue[ringStates, (ProjCOp[d, # + pauli] == -ProjCOp[d, #] && ProjSOp[d, # + pauli] == -ProjSOp[d, #]) &] &&
  AllTrue[ringStates, Module[{p = PythagoreanLiftOp[d, #]},
    p[[1]]^C0OverCSqOp[d] + p[[2]]^C0OverCSqOp[d] == p[[3]]^C0OverCSqOp[d]
  ] &] &&
  (Total[Table[ProjCOp[d, phi], {phi, ringStates}]] == 0) &&
  (Total[Table[ProjSOp[d, phi], {phi, ringStates}]] == 0) &&
  (Total[Table[ProjCOp[d, phi] * ProjSOp[d, phi], {phi, ringStates}]] == 0) &&
  (FresnelCosSumOp[d] == Quotient[KPOp[d]^C0OverCSqOp[d], C0OverCSqOp[d]]) &&
  (FresnelSinSumOp[d] == Quotient[KPOp[d]^C0OverCSqOp[d], C0OverCSqOp[d]])
];
