(* ::Package:: *)
(* Euclidean phase split: dphi == N_ring(d) * k_xy + A_xy; Maxwell B_f == N_ring(d) * n_f *)

BranchQuotientOp[d_Integer, dphi_Integer] :=
  Floor[(dphi + (PhiPauliOp[d] - MuVacOp[d])) / NRingOp[d]];

EdgeConnectionOp[d_Integer, dphi_Integer] :=
  dphi - NRingOp[d] * BranchQuotientOp[d, dphi];

EdgeDifferenceWithBranchCarryOp[d_Integer, phiY_Integer, phiX_Integer] := Module[
  {uY = Mod[phiY, NRingOp[d]], uX = Mod[phiX, NRingOp[d]], dphi},
  dphi = uY - uX;
  {EdgeConnectionOp[d, dphi], BranchQuotientOp[d, dphi]}
];

PlaquetteWindingOp[d_Integer, p00_Integer, p01_Integer, p11_Integer, p10_Integer] := Module[
  {k01, k12, k23, k30},
  k01 = EdgeDifferenceWithBranchCarryOp[d, p01, p00][[2]];
  k12 = EdgeDifferenceWithBranchCarryOp[d, p11, p01][[2]];
  k23 = EdgeDifferenceWithBranchCarryOp[d, p10, p11][[2]];
  k30 = EdgeDifferenceWithBranchCarryOp[d, p00, p10][[2]];
  -(k01 + k12 + k23 + k30)
];

PlaquetteMagneticFluxOp[d_Integer, p00_Integer, p01_Integer, p11_Integer, p10_Integer] := Module[
  {a01, a12, a23, a30},
  a01 = EdgeDifferenceWithBranchCarryOp[d, p01, p00][[1]];
  a12 = EdgeDifferenceWithBranchCarryOp[d, p11, p01][[1]];
  a23 = EdgeDifferenceWithBranchCarryOp[d, p10, p11][[1]];
  a30 = EdgeDifferenceWithBranchCarryOp[d, p00, p10][[1]];
  a01 + a12 + a23 + a30
];

(* Discrete Gauss monopole-free law over the 6 faces of a 3D unit cube: div_cube n_f === 0 *)
CubeMagneticDivergenceOp[d_Integer,
  p000_Integer, p001_Integer, p010_Integer, p011_Integer,
  p100_Integer, p101_Integer, p110_Integer, p111_Integer] := Module[
  {nZ0, nZ1, nY0, nY1, nX0, nX1},
  nZ0 = PlaquetteWindingOp[d, p000, p001, p011, p010];
  nZ1 = PlaquetteWindingOp[d, p100, p101, p111, p110];
  nY0 = PlaquetteWindingOp[d, p000, p001, p101, p100];
  nY1 = PlaquetteWindingOp[d, p010, p011, p111, p110];
  nX0 = PlaquetteWindingOp[d, p000, p010, p110, p100];
  nX1 = PlaquetteWindingOp[d, p001, p011, p111, p101];
  (nZ1 - nZ0) - (nY1 - nY0) + (nX1 - nX0)
];

VerifyHolonomyMaxwellModuleOp[d_Integer] := Module[
  {nRing = NRingOp[d], qSq = QSquareOp[d], mu = MuVacOp[d], pA, pB, pC, pD, w, flux},
  pA = 0;
  pB = qSq + C0OverCSqOp[d];
  pC = C0OverCSqOp[d] * pB;
  pD = d * pB;
  w = PlaquetteWindingOp[d, pA, pB, pC, pD];
  flux = PlaquetteMagneticFluxOp[d, pA, pB, pC, pD];
  (w == mu) &&
  (flux == nRing * w) &&
  (CubeMagneticDivergenceOp[d, pA, pB, pC, pD, pB, pC, pD, pA] == 0)
];
