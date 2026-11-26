(* ::Package:: *)
(* On-demand verification suite for Core and Lab packages *)

VerifySpecModuleOp[] := Module[{solutions, d},
  solutions = SolveUniqueCatalanDimensionOp[8];
  d = First[solutions];
  (Length[solutions] == MuVacOp[d]) &&
  CatalanLockPred[d] &&
  CatalanInductiveDominancePred[d + MuVacOp[d]] &&
  HeisenbergLockPred[d] &&
  GaussianCoulombLockPred[d] &&
  StarBooleanTwoLaneLockPred[d] &&
  (Length[AxisRotationPairsOp[d]] * C0OverCSqOp[d] == NNeighborsOp[d]) &&
  (Mod[DeltaPhiDiscOp[d] * BezoutInvDeltaPhiOp[d], NRingOp[d]] == MuVacOp[d]) &&
  (ShellSecondMomentTensorOp[d, O2CuboctahedronVerticesOp[d]] ==
    NTrianglesOp[d] * IdentityMatrix[d]) &&
  (ShellCrossFourthMomentOp[O1SquareFaceVerticesOp[d]] == 0) &&
  (ShellCrossFourthMomentOp[O2CuboctahedronVerticesOp[d]] == DInSquareOp[d]) &&
  (!ShellHasPlanarZeroCoordPred[O3CornerOctantVerticesOp[d]]) &&
  ShellHasPlanarZeroCoordPred[O2CuboctahedronVerticesOp[d]]
];

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

VerifyCarrierAutomaton3DModuleOp[d_Integer] := Module[
  {two, g0, gWallPos, gWallNeg, ivGrid, inertGrid, inert1, bzGrid0, bzGrid1, bzRev, cIdx},
  two      = C0OverCSqOp[d];
  g0       = Compute12NeighborHolonomyGateAtIvOp[d, two, two, two, 0];
  gWallPos = Wrap512Op[d, Compute12NeighborHolonomyGateAtIvOp[d, 0, two, two, MuVacOp[d]]];
  gWallNeg = Wrap512Op[d, Compute12NeighborHolonomyGateAtIvOp[d, -MuVacOp[d], two, two, MuVacOp[d]]];

  ivGrid    = CreateIvCarrier3DOp[d, d + MuVacOp[d]];
  inertGrid = CreateReactionCarrier3DOp[d, d, DeltaPhiDiscOp[d]];
  inert1    = StepCarrierForward3DOp[d, inertGrid];
  cIdx      = two;

  bzGrid0 = CreateReactionCarrier3DOp[d, d, CollapsePhiOp[d]];
  bzGrid1 = StepCarrierForward3DOp[d, bzGrid0];
  bzRev   = StepCarrierBackward3DOp[d, bzGrid1];

  (RingVisibilityL1NormOp[d, 0, DeltaPhiDiscOp[d]] == 0) &&
  (RingVisibilityL1NormOp[d, MuVacOp[d], DeltaPhiDiscOp[d]] == QSquareOp[d]) &&
  (EnforcePlanckCellFloorOp[d, 0] == MuVacOp[d]) &&
  AllCellsFilledPred[d, ivGrid] &&
  AllCellsFilledPred[d, bzGrid1] &&
  (g0 == 0) &&
  (gWallPos + gWallNeg == -MuVacOp[d]) &&
  (inert1[[cIdx, cIdx, cIdx, 2]] == 0) &&
  (bzGrid1[[cIdx, cIdx, cIdx, 2]] == NRingOp[d] - DOutEdgeOp[d] * Quotient[QSquareOp[d], DInSquareOp[d]]) &&
  (bzRev[[All, All, All, 1 ;; 6]] == bzGrid0[[All, All, All, 1 ;; 6]])
];

VerifyContinuumPhysicsModuleOp[d_Integer] := (
  (EnergyQuantumOp[d, PhiPauliOp[d]] == NSquaresOp[d]) &&
  (EnergyQuantumOp[d, DeltaPhiDiscOp[d] - MuVacOp[d]] == 0) &&
  (AlphaFloorIntOp[d] == (NPhiOp[d] - C0OverCSqOp[d])^C0OverCSqOp[d] + DInSquareOp[d]^C0OverCSqOp[d]) &&
  (AlphaFloorIntOp[d] == QSquareOp[d] + ModBitsOp[d]) &&
  (AlphaFloorLhsOp[d, 0]^C0OverCSqOp[d] < AlphaFloorRhsSqOp[d]) &&
  (AlphaFloorRhsSqOp[d] < AlphaFloorLhsOp[d, MuVacOp[d]]^C0OverCSqOp[d]) &&
  (ShellGateDenomOp[d, KPOp[d], 0] == MaxGravitationalCompressionRatioOp[d] * EpsDenomOp[d]) &&
  (ShellGateDenomOp[d, KPOp[d], CollapsePhiOp[d]] == EpsDenomOp[d]) &&
  (ShellGateDenomOp[d, KPOp[d], PhiPauliOp[d]] == EpsDenomOp[d]) &&
  (DarkMatterResidualOp[d, NPhiOp[d], DeltaPhiDiscOp[d]] == 0) &&
  (VoronoiHaloToCoreRatioOp[d] == (C0OverCSqOp[d] * NTrianglesOp[d]) / d) &&
  (TsirelsonQmSqOp[d] == NTrianglesOp[d])
);

VerifyBioCosmologyModuleOp[d_Integer] := Module[
  {cageN4 = TetrapyrroleCageChargeOp[d], cageN6 = HexahistidineCageChargeOp[d],
   hb = HBondDivisorOp[d], mu = MuVacOp[d], two = C0OverCSqOp[d]},
  (CodonSpaceOp[d] == KPOp[d]) &&
  (two * AminoAcidClassesOp[d] == CuboctShellSizeOp[d, two]) &&
  (two * hb == d * CuboctShellSizeOp[d, two]) &&
  (CuboctBallOp[d] == d * (NTrianglesOp[d] - mu)^two) &&
  (cageN4 == NuclearMagicHighOp[d, DInSquareOp[d]]) &&
  (cageN6 == CuboctShellSizeOp[d, two]) &&
  (DiscoverBioMetalInCageOp[KPOp[d], cageN4, 0, DeltaPhiDiscOp[d] - mu] == {ZMagnesiumOp[d]}) &&
  (DiscoverBioMetalInCageOp[KPOp[d], cageN4, 0, CuboctBall2Op[d]] == {d^d}) &&
  (DiscoverBioMetalInCageOp[KPOp[d], cageN4, 0, two * NuclearMagicHighOp[d, DInSquareOp[d]]] == {ZNickelOp[d]}) &&
  (DiscoverBioMetalInCageOp[KPOp[d], cageN4, two * ZOxygenOp[d], KPOp[d] + NSquaresOp[d]] == {ZIronOp[d]}) &&
  (DiscoverBioMetalInCageOp[KPOp[d], ZSulfurOp[d], 0, CuboctShellSizeOp[d, two]] == {ZIronOp[d]}) &&
  (DiscoverBioMetalInCageOp[KPOp[d], two * ZSulfurOp[d], two * ZNitrogenOp[d], KPOp[d] + NNeighborsOp[d]] == {d * ZNeonOp[d]}) &&
  ({UniversalPhotonStepOp[d, -mu] / hb, UniversalPhotonStepOp[d, 0] / hb, UniversalPhotonStepOp[d, mu] / hb} ==
   {ZNitrogenOp[d], two, mu})
];

VerifyNuclearShellsModuleOp[d_Integer] := Module[
  {nCoding = AminoAcidClassesOp[d] - MuVacOp[d], codingRange, two = C0OverCSqOp[d], dphi = DeltaPhiDiscOp[d]},
  codingRange = Range[0, nCoding - 1];
  (Period8RemainingElementsOp[d] == NPhiOp[d] + NSquaresOp[d]) &&
  (NIslandOfStabilityOp[d] == two * ZUraniumOp[d]) &&
  (ZLeadTerminalOp[d] < ZUraniumOp[d] < ZIslandFlOp[d] < ZOganessonOp[d] <
   ZIslandUbnOp[d] < ZCoreDiracDivingOp[d] < ZFissionCeilingOp[d] <
   ZFeynmanPointLimitOp[d] < ZQedFiniteNucleusLowOp[d] < ZQedFiniteNucleusHighOp[d] <
   NIslandOfStabilityOp[d]) &&
  AllTrue[codingRange, (TripletCodonToAminoAcidOp[d, AminoAcidToTripletCodonOp[d, #]] == #) &] &&
  AllTrue[Tuples[{codingRange, codingRange}],
    IsJunctionFrameProtectedPred[d, AminoAcidToTripletCodonOp[d, #[[1]]], AminoAcidToTripletCodonOp[d, #[[2]]]] &
  ] &&
  (TissueEuler2DOp[GrowTissueFromDnaSequenceInEnvOp[d, CompilePhenotypeToDnaSequenceOp[d, two], 0, 0, dphi, d]] == two) &&
  (TissueEuler2DOp[GrowTissueFromDnaSequenceInEnvOp[d, CompilePhenotypeToDnaSequenceOp[d, 0], 0, 0, dphi, d]] == 0) &&
  (TissueEuler2DOp[GrowTissueFromDnaSequenceInEnvOp[d, CompilePhenotypeToDnaSequenceOp[d, dphi], 0, 0, dphi, d]] == dphi) &&
  (DnaStringToIndicesOp[d, "ATG"] == {two, 0, d}) &&
  (TissueEuler2DOp[GrowTissueFromDnaSequenceInEnvOp[d, CompilePhenotypeToDnaStringOp[d, dphi], 0, 0, dphi, d]] == dphi) &&
  (TissueEuler2DOp[GrowTissueFromDnaSequenceInEnvOp[d, CompilePhenotypeToBioSequenceOp[d, dphi], 0, 0, dphi, d]] == dphi)
];

VerifyParticleMaterialsModuleOp[d_Integer] := Module[
  {mu = MuVacOp[d], two = C0OverCSqOp[d], qSq = QSquareOp[d], pauli = PhiPauliOp[d],
   dphiDisc = DeltaPhiDiscOp[d], nPhi = NPhiOp[d], aCoh = CollapsePhiOp[d]},
  (BondEnergyOp[d, mu, mu, mu, dphiDisc] + BondBohmDeficitOp[d, mu, mu, mu, dphiDisc] == qSq) &&
  (ProjCOp[d, qSq] == 0 && Abs[ProjSOp[d, qSq]] == qSq) &&
  (ProjCOp[d, d * dphiDisc] == DOutEdgeOp[d]) &&
  (ProjCOp[d, DInSquareOp[d] * dphiDisc] == -SingleTickSpinorBitBudgetOp[d]) &&
  (CatalyticCavityGateOp[d, nPhi, dphiDisc] == 0) &&
  (Abs[CatalyticCavityGateOp[d, nPhi, aCoh]] > qSq) &&
  (CooperPairProjCOp[d, qSq] == -qSq && CooperPairProjSOp[d, qSq] == 0) &&
  (CooperPairProjCOp[d, pauli] == qSq && CooperPairProjSOp[d, pauli] == 0) &&
  (ZNitinolShapeMemoryOp[d] == NuclearMagicHighOp[d, DOutEdgeOp[d]]) &&
  (ZMagnesiumDiborideOp[d] == ZTitaniumBeamOp[d]) &&
  (ZCupratePlaneOp[d] == ModBitsOp[d] * DOutEdgeOp[d]) &&
  (ZInfiniteLayerCuprateOp[d] == KPOp[d] + mu) &&
  (ZIronSulfideSuperOp[d] == CuboctShellSizeOp[d, two])
];

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

VerifyAllCarrierModulesOp[d_Integer] := {
  {"1. Core/Spec.wl",                  VerifySpecModuleOp[]},
  {"2. Core/FaceAlgebra.wl",           VerifyFaceAlgebraModuleOp[d]},
  {"3. Core/HolonomyMaxwell.wl",       VerifyHolonomyMaxwellModuleOp[d]},
  {"4. Core/CarrierAutomaton3D.wl",    VerifyCarrierAutomaton3DModuleOp[d]},
  {"5. Lab/ContinuumPhysics.wl",       VerifyContinuumPhysicsModuleOp[d]},
  {"6. Lab/BioCosmology.wl",           VerifyBioCosmologyModuleOp[d]},
  {"7. Lab/NuclearShells.wl",          VerifyNuclearShellsModuleOp[d]},
  {"8. Lab/ParticleMaterials.wl",      VerifyParticleMaterialsModuleOp[d]},
  {"9. Lab/ReversibleCalculus.wl",     VerifyReversibleCalculusModuleOp[d]},
  {"10. Lab/HardwarePhotovoltaics.wl", VerifyHardwarePhotovoltaicsModuleOp[d]}
};
