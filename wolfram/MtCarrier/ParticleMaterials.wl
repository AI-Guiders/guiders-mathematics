(* ::Package:: *)
(* Covalent bond holonomy, saddle |A^*| == Q_sq(d), catalytic gate on O_2(d), Cooper 2e *)

BondEnergyOp[d_Integer, b_Integer, muX_Integer, muY_Integer, dphi_Integer] :=
  b * (muX * muY) * ProjCOp[d, dphi];

BondBohmDeficitOp[d_Integer, b_Integer, muX_Integer, muY_Integer, dphi_Integer] :=
  b * (muX * muY) * Abs[Wrap512Op[d, dphi]];

ReactionCurrentOp[d_Integer, b_Integer, muX_Integer, muY_Integer, dphi_Integer] :=
  b * (muX * muY) * ProjSOp[d, dphi];

ReactionEnthalpyOp[hReactants_Integer, hProducts_Integer] :=
  hReactants - hProducts;

CatalyticCavityGateOp[d_Integer, muReactant_Integer, dphiEnv_Integer] := Module[
  {sumMu = NNeighborsOp[d] * muReactant, zr, zi = 0},
  zr = NNeighborsOp[d] * (muReactant^C0OverCSqOp[d]) * ProjCOp[d, dphiEnv];
  GatePhaseOp[d, muReactant, sumMu, zr, zi]
];

CooperPairProjCOp[d_Integer, dphi_Integer] :=
  ProjCOp[d, C0OverCSqOp[d] * dphi];

CooperPairProjSOp[d_Integer, dphi_Integer] :=
  ProjSOp[d, C0OverCSqOp[d] * dphi];

ZTitaniumBeamOp[d_Integer] := C0OverCSqOp[d] * (NNeighborsOp[d] - MuVacOp[d]);
ZCopperOp[d_Integer]       := DInSquareOp[d] * DOutVertexOp[d] + MuVacOp[d];

ZNitinolShapeMemoryOp[d_Integer]   := ZTitaniumBeamOp[d] + ZNickelOp[d];
ZMagnesiumDiborideOp[d_Integer]    := ZMagnesiumOp[d] + C0OverCSqOp[d] * ZBoronOp[d];
ZCupratePlaneOp[d_Integer]         := ZCopperOp[d] + C0OverCSqOp[d] * ZOxygenOp[d];
ZInfiniteLayerCuprateOp[d_Integer] := ZCalciumOp[d] + ZCupratePlaneOp[d];
ZIronSulfideSuperOp[d_Integer]     := ZIronOp[d] + ZSulfurOp[d];

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
