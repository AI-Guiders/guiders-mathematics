(* ::Package:: *)
(* Lab Layer above Core CA: Algebraic alpha_fs^-1(d), gravity floor eps(d), Voronoi halo *)

AlphaNumAOp[d_Integer] :=
  Quotient[DInSquareOp[d], C0OverCSqOp[d]] * (NTrianglesOp[d] - MuVacOp[d]) * (MHatOp[d]^C0OverCSqOp[d]);

AlphaNumBOp[d_Integer] :=
  Quotient[DInSquareOp[d], C0OverCSqOp[d]] * ((NTrianglesOp[d] - MuVacOp[d])^C0OverCSqOp[d]) * (MHatOp[d]^d);

AlphaDenAOp[d_Integer] :=
  DInSquareOp[d] * (((NTrianglesOp[d] - MuVacOp[d])^C0OverCSqOp[d]) * (MHatOp[d]^DInSquareOp[d]) - NTrianglesOp[d]);

AlphaDenBOp[d_Integer] :=
  -(NTrianglesOp[d] - MuVacOp[d]) * MHatOp[d];

AlphaRatDenOp[d_Integer] :=
  C0OverCSqOp[d] * (AlphaNumBOp[d]^C0OverCSqOp[d]) - (AlphaNumAOp[d]^C0OverCSqOp[d]);

AlphaRatPOp[d_Integer] :=
  C0OverCSqOp[d] * AlphaDenBOp[d] * AlphaNumBOp[d] - AlphaDenAOp[d] * AlphaNumAOp[d];

AlphaRatQOp[d_Integer] :=
  AlphaDenAOp[d] * AlphaNumBOp[d] - AlphaDenBOp[d] * AlphaNumAOp[d];

AlphaFloorIntOp[d_Integer] :=
  MHatOp[d] + DOutEdgeOp[d] * NTrianglesOp[d];

AlphaFloorLhsOp[d_Integer, k_Integer] :=
  (AlphaFloorIntOp[d] + k) * AlphaRatDenOp[d] - AlphaRatPOp[d];

AlphaFloorRhsSqOp[d_Integer] :=
  C0OverCSqOp[d] * (AlphaRatQOp[d]^C0OverCSqOp[d]);

AlphaInvExactAlgebraicOp[d_Integer] :=
  (AlphaRatPOp[d] + AlphaRatQOp[d] * Sqrt[C0OverCSqOp[d]]) / AlphaRatDenOp[d];

ShellGateDenomOp[d_Integer, mu_Integer, dphi_Integer] :=
  GateDenomOp[d, mu, NNeighborsOp[d] * (mu^C0OverCSqOp[d]) * ProjCOp[d, dphi]];

MaxGravitationalCompressionRatioOp[d_Integer] :=
  MuVacOp[d] + NPhiOp[d] * KPOp[d];

BaryonStaticDenomOp[d_Integer, mu_Integer] :=
  NPhiOp[d] * QSquareOp[d] * (mu^C0OverCSqOp[d]) + EpsDenomOp[d];

BohmGaugeDeficitOp[d_Integer, mu_Integer, dphi_Integer] :=
  NNeighborsOp[d] * (mu^C0OverCSqOp[d]) * Abs[Wrap512Op[d, dphi]];

DarkMatterResidualOp[d_Integer, mu_Integer, dphi_Integer] :=
  (NNeighborsOp[d] * (mu^C0OverCSqOp[d]) * ProjCOp[d, dphi] + QSquareOp[d] * (mu^C0OverCSqOp[d]) + EpsDenomOp[d]) -
  (BaryonStaticDenomOp[d, mu] - BohmGaugeDeficitOp[d, mu, dphi]);

VoronoiHaloToCoreRatioOp[d_Integer] :=
  (C0OverCSqOp[d]^ModBitsOp[d] * DInSquareOp[d]) / (d * QSquareOp[d]);

BellProdSqOp[d_Integer]    := C0OverCSqOp[d];
BellLhvSqOp[d_Integer]     := DInSquareOp[d];
TsirelsonQmSqOp[d_Integer] := C0OverCSqOp[d] * BellLhvSqOp[d];
