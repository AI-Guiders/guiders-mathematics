(* ::Package:: *)
(* Lab Layer above Core CA: 12-beam splitter, 4-channel Neumann tomograph, singlet fission *)

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

BluePhotonWavelengthNmOp[d_Integer]    := SeamTicksOp[d]^C0OverCSqOp[d];
ExcitonPairWavelengthNmOp[d_Integer]   := C0OverCSqOp[d] * BluePhotonWavelengthNmOp[d];
SiliconBandgapCutoffNmOp[d_Integer]    := (d^d) * DeltaPhiDiscOp[d];
TetraceneCoreChargeOp[d_Integer]       := NNeighborsOp[d] * ZNeonOp[d];
RubreneQuinoneLockChargeOp[d_Integer]  := TetraceneCoreChargeOp[d] + ZOxygenOp[d];
