(* ::Package:: *)
(* 12-beam splitter, 4-channel Neumann tomograph (2 I_4 + N_nil)^-1, silicon singlet fission *)

TLayerHiddenMicrostatesPerCellOp[d_Integer] :=
  KPOp[d] * HeisenbergEquivClassSizeOp[d];

MaxwellDemonRayPhaseOp[d_Integer, phiThermal_Integer] :=
  Quotient[phiThermal, NNeighborsOp[d]];

MaxwellDemonTotalExtractedCurrentOp[d_Integer, phiThermal_Integer] :=
  NNeighborsOp[d] * ProjSOp[d, MaxwellDemonRayPhaseOp[d, phiThermal]];

NilpotentTomographForward2xOp[d_Integer, {s0_Integer, s1_Integer, s2_Integer, s3_Integer}] := Module[
  {two = C0OverCSqOp[d]},
  {two * s0 + s1, two * s1 + two * s2, two * s2, two * s3}
];

NilpotentTomographInverse24xOp[d_Integer, {y0_Integer, y1_Integer, y2_Integer, y3_Integer}] := Module[
  {n12 = NNeighborsOp[d], n6 = NSquaresOp[d], den = C0OverCSqOp[d] * NNeighborsOp[d]},
  {Quotient[n12 * y0 - n6 * y1 + n6 * y2, den],
   Quotient[n12 * y1 - n12 * y2, den],
   Quotient[n12 * y2, den],
   Quotient[n12 * y3, den]}
];

(* Singlet fission: Seam(d)^2 nm -> 2 x (2 * Seam(d)^2) nm < d^d * DPhiDisc(d) nm *)
BluePhotonWavelengthNmOp[d_Integer]    := SeamTicksOp[d]^C0OverCSqOp[d];
ExcitonPairWavelengthNmOp[d_Integer]   := C0OverCSqOp[d] * BluePhotonWavelengthNmOp[d];
SiliconBandgapCutoffNmOp[d_Integer]    := (d^d) * DeltaPhiDiscOp[d];
TetraceneCoreChargeOp[d_Integer]       := NNeighborsOp[d] * ZNeonOp[d];
RubreneQuinoneLockChargeOp[d_Integer]  := TetraceneCoreChargeOp[d] + ZOxygenOp[d];

VerifyHardwarePhotovoltaicsModuleOp[d_Integer] := Module[
  {pauli = PhiPauliOp[d], sampleVec, mixedVec, recVec},
  sampleVec = {ZSiliconOp[d], ZSulfurOp[d], d * ZNeonOp[d], ZSiliconOp[d] + DInSquareOp[d] * ZOxygenOp[d]};
  mixedVec  = NilpotentTomographForward2xOp[d, sampleVec];
  recVec    = NilpotentTomographInverse24xOp[d, mixedVec];

  (TLayerHiddenMicrostatesPerCellOp[d] == (TwoTickPhaseSpaceBitBudgetOp[d]^C0OverCSqOp[d])) &&
  (MaxwellDemonRayPhaseOp[d, pauli] == SeamTicksOp[d]) &&
  (EnergyQuantumOp[d, MaxwellDemonRayPhaseOp[d, pauli]] == 0) &&
  (MaxwellDemonTotalExtractedCurrentOp[d, pauli] == pauli - DInSquareOp[d]) &&
  (recVec == sampleVec) &&
  (BluePhotonWavelengthNmOp[d] < ExcitonPairWavelengthNmOp[d] < SiliconBandgapCutoffNmOp[d]) &&
  (RubreneQuinoneLockChargeOp[d] == QSquareOp[d])
];
