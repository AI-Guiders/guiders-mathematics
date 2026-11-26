(* ::Package:: *)
(* Zero-Landauer-heat Verlet, |E_Heis|(d) == d^(d_in(sq)), Born rule P(mu) == mu^2 / N_ring(d), 3-SAT *)

HeisenbergEquivClassSizeOp[d_Integer] :=
  C0OverCSqOp[d] * DeltaPhiDiscOp[d] - MuVacOp[d];

DiracBraKetOp[d_Integer, mu1_Integer, phi1_Integer, mu2_Integer, phi2_Integer] :=
  {mu1 * mu2 * ProjCOp[d, phi2 - phi1],
   mu1 * mu2 * ProjSOp[d, phi2 - phi1]};

(* Counts background detector phases phiBg in Z_{N_ring(d)} where branch quotient click == mu_vac(d) *)
BornClickCountOp[d_Integer, mu_Integer] := Module[
  {phiWave = mu^C0OverCSqOp[d], muVac = MuVacOp[d]},
  Sum[
    EqIntIndicatorOp[BranchQuotientOp[d, Wrap512Op[d, phiBg] + phiWave], muVac],
    {phiBg, 0, NRingOp[d] - muVac}
  ]
];

Clause3SatDefectOp[d_Integer, psi1_Integer, psi2_Integer, psi3_Integer] :=
  QSquareOp[d] - (ProjCOp[d, psi1] + ProjCOp[d, psi2] + ProjCOp[d, psi3]);

Clause3SatVortexQuantumOp[d_Integer, psi1_Integer, psi2_Integer, psi3_Integer] :=
  IntegerPart[Clause3SatDefectOp[d, psi1, psi2, psi3] / NRingOp[d]];

Clause3AndQuantumOp[d_Integer, psi1_Integer, psi2_Integer, psi3_Integer] :=
  IntegerPart[(QSquareOp[d] + (ProjCOp[d, psi1] + ProjCOp[d, psi2] + ProjCOp[d, psi3])) / NRingOp[d]];

Count3SatZeroQuotientOctantsOp[d_Integer] := Module[
  {poles = {0, PhiPauliOp[d]}, octants},
  octants = Tuples[ConstantArray[poles, d]];
  Sum[
    IsZeroNatOp[Clause3SatVortexQuantumOp[d, octants[[k, 1]], octants[[k, 2]], octants[[k, 3]]]],
    {k, 1, Length[octants]}
  ]
];

CuboctCausalConeVolumeOp[d_Integer, tDepth_Integer] :=
  Sum[CuboctShellSizeOp[d, k], {k, 0, tDepth}];

LocalInversionPreimageBranchingOp[d_Integer, tDepth_Integer] :=
  (NTrianglesOp[d] - MuVacOp[d]) ^ tDepth;

VerifyReversibleCalculusModuleOp[d_Integer] := Module[
  {mu, two, maxBornMu},
  mu = MuVacOp[d];
  two = C0OverCSqOp[d];
  maxBornMu = two * NTrianglesOp[d];
  (HeisenbergEquivClassSizeOp[d] == (d ^ DInSquareOp[d])) &&
  (HeisenbergEquivClassSizeOp[d] == (ModBitsOp[d] ^ two)) &&
  (CuboctShellSizeOp[d, DInSquareOp[d]] == two * HeisenbergEquivClassSizeOp[d]) &&
  AllTrue[Range[0, maxBornMu], (BornClickCountOp[d, #] == (#^two)) &] &&
  (Count3SatZeroQuotientOctantsOp[d] == NTrianglesOp[d] - mu) &&
  (CuboctCausalConeVolumeOp[d, d] == CuboctBallOp[d]) &&
  (LocalInversionPreimageBranchingOp[d, d] == (DOutVertexOp[d] ^ d))
];
