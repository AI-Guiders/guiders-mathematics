(* ::Package:: *)
(* Pure M-Layer Core: Discrete Euler character C_d, S_d on Z_{N_ring(d)} *)

IsZeroNatOp[n_Integer] := Quotient[1, Abs[n] + 1];
IsPosNatOp[n_Integer]  := 1 - IsZeroNatOp[Max[0, n]];
EqIntIndicatorOp[a_Integer, b_Integer] := IsZeroNatOp[Abs[a - b]];
BlendByIndicatorOp[sel_Integer, valTrue_, valFalse_] := sel * valTrue + (1 - sel) * valFalse;

Wrap512Op[d_Integer, phi_Integer] := Module[
  {u = Mod[phi, NRingOp[d]], q},
  q = Quotient[u + (PhiPauliOp[d] - MuVacOp[d]), NRingOp[d]];
  u - NRingOp[d] * q
];

RingFaceLutCOp[d_Integer] := RingFaceLutCOp[d] =
  Table[QSquareOp[d] - Abs[Wrap512Op[d, p]], {p, 0, NRingOp[d] - MuVacOp[d]}];

RingFaceLutSOp[d_Integer] := RingFaceLutSOp[d] =
  Table[RingFaceLutCOp[d][[Mod[p - QSquareOp[d], NRingOp[d]] + MuVacOp[d]]], {p, 0, NRingOp[d] - MuVacOp[d]}];

ProjCOp[d_Integer, phi_Integer] := RingFaceLutCOp[d][[Mod[phi, NRingOp[d]] + MuVacOp[d]]];
ProjSOp[d_Integer, phi_Integer] := RingFaceLutSOp[d][[Mod[phi, NRingOp[d]] + MuVacOp[d]]];

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
