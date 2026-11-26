(* ::Package:: *)
(* Root entry point mirroring lean/MtCarrier.lean: loads all 10 MtCarrier/*.wl modules *)

mtCarrierBaseDir = DirectoryName[$InputFileName];
mtCarrierPkgDir  = FileNameJoin[{mtCarrierBaseDir, "MtCarrier"}];

Get[FileNameJoin[{mtCarrierPkgDir, "Spec.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "FaceAlgebra.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "HolonomyMaxwell.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "CarrierAutomaton3D.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "ContinuumPhysics.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "BioCosmology.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "NuclearShells.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "ParticleMaterials.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "ReversibleCalculus.wl"}]];
Get[FileNameJoin[{mtCarrierPkgDir, "HardwarePhotovoltaics.wl"}]];

dSolutions = SolveUniqueCatalanDimensionOp[8];
dSolved    = First[dSolutions];

Print["======================================================================"];
Print["MODULAR WOLFRAM LANGUAGE DQM SUITE (10 PACKAGES, 0 HARDCODED CONSTANTS)"];
Print["======================================================================"];
Print["Solved unique Catalan dimension d from |CellBranchingStates(d)| == |AxisTransitionPairs(d)|:"];
Print["  Solutions in 1..8: ", dSolutions, " -> d = ", dSolved];
Print["  Derived Passport from d alone:"];
Print["    mu_vac = ", MuVacOp[dSolved],
      ", (c0/c)^2 = ", C0OverCSqOp[dSolved],
      ", |O_2(d)| = ", NNeighborsOp[dSolved],
      ", |tri(d)| = ", NTrianglesOp[dSolved],
      ", |sq(d)| = ", NSquaresOp[dSolved]];
Print["    mod_bits = ", ModBitsOp[dSolved],
      ", N_ring = ", NRingOp[dSolved],
      ", K_P = ", KPOp[dSolved],
      ", Q_sq = ", QSquareOp[dSolved],
      ", Phi_Pauli = ", PhiPauliOp[dSolved]];
Print["    N_phi = ", NPhiOp[dSolved],
      ", DeltaPhi_disc = ", DeltaPhiDiscOp[dSolved],
      ", Seam = ", SeamTicksOp[dSolved],
      ", BezoutInv = ", BezoutInvDeltaPhiOp[dSolved],
      ", M_hat = ", MHatOp[dSolved],
      ", IV(d) = ", IvSeedOp[dSolved]];

suiteResults = {
  {"1. MtCarrier/Spec.wl",                  VerifySpecModuleOp[]},
  {"2. MtCarrier/FaceAlgebra.wl",           VerifyFaceAlgebraModuleOp[dSolved]},
  {"3. MtCarrier/HolonomyMaxwell.wl",       VerifyHolonomyMaxwellModuleOp[dSolved]},
  {"4. MtCarrier/CarrierAutomaton3D.wl",    VerifyCarrierAutomaton3DModuleOp[dSolved]},
  {"5. MtCarrier/ContinuumPhysics.wl",      VerifyContinuumPhysicsModuleOp[dSolved]},
  {"6. MtCarrier/BioCosmology.wl",          VerifyBioCosmologyModuleOp[dSolved]},
  {"7. MtCarrier/NuclearShells.wl",         VerifyNuclearShellsModuleOp[dSolved]},
  {"8. MtCarrier/ParticleMaterials.wl",     VerifyParticleMaterialsModuleOp[dSolved]},
  {"9. MtCarrier/ReversibleCalculus.wl",    VerifyReversibleCalculusModuleOp[dSolved]},
  {"10. MtCarrier/HardwarePhotovoltaics.wl", VerifyHardwarePhotovoltaicsModuleOp[dSolved]}
};

Do[
  Print["  ", suiteResults[[k, 1]], " -> ", suiteResults[[k, 2]]],
  {k, 1, Length[suiteResults]}
];

allModulesOk = AllTrue[suiteResults, #[[2]] === True &];
Print["======================================================================"];
Print["ALL 10 MODULAR WOLFRAM LANGUAGE PACKAGES VERIFIED: ", allModulesOk];
Print["======================================================================"];
