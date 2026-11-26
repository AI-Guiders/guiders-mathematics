(* ::Package:: *)
(* Pure M-Layer Core: Catalan closure 1 + 2^d == d^2 -> unique d -> O_2(d) *)

CentralOccupancyStateOp[d_Integer] := CentralOccupancyStateOp[d] = ConstantArray[0, d];
ElementaryAxisSignsOp[d_Integer]   := ElementaryAxisSignsOp[d]   = {-Quotient[d, d], Quotient[d, d]};
SignedOctantsOp[d_Integer]         := SignedOctantsOp[d]         = Tuples[ConstantArray[ElementaryAxisSignsOp[d], d]];
CellBranchingStatesOp[d_Integer]   := CellBranchingStatesOp[d]   = Prepend[SignedOctantsOp[d], CentralOccupancyStateOp[d]];

AxisTransitionPairsOp[d_Integer]   := AxisTransitionPairsOp[d]   = Tuples[{Range[d], Range[d]}];
AxisLongitudinalPairsOp[d_Integer] := AxisLongitudinalPairsOp[d] = Select[AxisTransitionPairsOp[d], #[[1]] == #[[2]] &];
AxisRotationPairsOp[d_Integer]     := AxisRotationPairsOp[d]     = Select[AxisTransitionPairsOp[d], #[[1]] != #[[2]] &];

CatalanLockPred[d_Integer] :=
  Length[CellBranchingStatesOp[d]] == Length[AxisTransitionPairsOp[d]];

SolveUniqueCatalanDimensionOp[maxSearch_Integer] :=
  Select[Range[1, maxSearch], CatalanLockPred];

CatalanInductiveDominancePred[d_Integer] :=
  (Length[SignedOctantsOp[d]] >= Length[AxisTransitionPairsOp[d]]) &&
  (Length[SignedOctantsOp[d + 1]] > Length[AxisTransitionPairsOp[d + 1]]);

MuVacOp[d_Integer]    := MuVacOp[d]    = Length[AxisTransitionPairsOp[d]] - Length[SignedOctantsOp[d]];
C0OverCSqOp[d_Integer] := C0OverCSqOp[d] = d - MuVacOp[d];

MooreDisplacementCubeOp[d_Integer] := MooreDisplacementCubeOp[d] =
  Select[Tuples[ConstantArray[{-MuVacOp[d], 0, MuVacOp[d]}, d]], (# . #) > 0 &];

CubicOrbitByNormSqOp[d_Integer, normSq_Integer] := CubicOrbitByNormSqOp[d, normSq] =
  Select[MooreDisplacementCubeOp[d], (# . #) == normSq &];

O1SquareFaceVerticesOp[d_Integer]    := O1SquareFaceVerticesOp[d]    = CubicOrbitByNormSqOp[d, MuVacOp[d]];
O2CuboctahedronVerticesOp[d_Integer] := O2CuboctahedronVerticesOp[d] = CubicOrbitByNormSqOp[d, C0OverCSqOp[d]];
O3CornerOctantVerticesOp[d_Integer]  := O3CornerOctantVerticesOp[d]  = CubicOrbitByNormSqOp[d, d];

NNeighborsOp[d_Integer] := NNeighborsOp[d] = Length[O2CuboctahedronVerticesOp[d]];
NTrianglesOp[d_Integer] := NTrianglesOp[d] = Length[SignedOctantsOp[d]];
NSquaresOp[d_Integer]   := NSquaresOp[d]   = Length[O1SquareFaceVerticesOp[d]];
DInSquareOp[d_Integer]  := DInSquareOp[d]  = Quotient[NNeighborsOp[d], d];
DOutEdgeOp[d_Integer]   := DOutEdgeOp[d]   = d + C0OverCSqOp[d];
DOutVertexOp[d_Integer] := DOutVertexOp[d] = NSquaresOp[d] + MuVacOp[d];

ModBitsOp[d_Integer]    := ModBitsOp[d]    = Length[CellBranchingStatesOp[d]];
FracBitsOp[d_Integer]   := FracBitsOp[d]   = ModBitsOp[d] - d;
NRingOp[d_Integer]      := NRingOp[d]      = C0OverCSqOp[d] ^ ModBitsOp[d];
KPOp[d_Integer]         := KPOp[d]         = C0OverCSqOp[d] ^ FracBitsOp[d];
QSquareOp[d_Integer]    := QSquareOp[d]    = Quotient[NRingOp[d], DInSquareOp[d]];
PhiPauliOp[d_Integer]   := PhiPauliOp[d]   = Quotient[NRingOp[d], C0OverCSqOp[d]];
NPhiOp[d_Integer]       := NPhiOp[d]       = MuVacOp[d] + NNeighborsOp[d];

HeisenbergDividendOp[d_Integer] := HeisenbergDividendOp[d] = NRingOp[d] + NPhiOp[d] + NTrianglesOp[d];
HeisenbergLockPred[d_Integer]   := Mod[HeisenbergDividendOp[d], NPhiOp[d]] == 0;
DeltaPhiDiscOp[d_Integer]       := DeltaPhiDiscOp[d]       = Quotient[HeisenbergDividendOp[d], NPhiOp[d]];
SeamTicksOp[d_Integer]          := SeamTicksOp[d]          = NPhiOp[d] * DeltaPhiDiscOp[d] - NRingOp[d];
BezoutInvDeltaPhiOp[d_Integer]  := BezoutInvDeltaPhiOp[d]  = PowerMod[DeltaPhiDiscOp[d], -MuVacOp[d], NRingOp[d]];

MHatOp[d_Integer] := MHatOp[d] = MuVacOp[d] + NNeighborsOp[d] * NTrianglesOp[d];
GaussianCoulombLockPred[d_Integer] :=
  MHatOp[d] == (DInSquareOp[d]^C0OverCSqOp[d] + ModBitsOp[d]^C0OverCSqOp[d]);

StarBooleanTwoLaneLockPred[d_Integer] :=
  (C0OverCSqOp[d] * NSquaresOp[d] == NNeighborsOp[d]) &&
  (C0OverCSqOp[d]^NNeighborsOp[d] == KPOp[d]^C0OverCSqOp[d]);

CollapsePhiOp[d_Integer] := CollapsePhiOp[d] = QSquareOp[d] + NNeighborsOp[d] - MuVacOp[d];
IvSeedOp[d_Integer]      := IvSeedOp[d]      = PhiPauliOp[d] + (NTrianglesOp[d] - MuVacOp[d]);
EpsDenomOp[d_Integer]    := EpsDenomOp[d]    = KPOp[d] * QSquareOp[d];

SingleTickSpinorBitBudgetOp[d_Integer]  := SingleTickSpinorBitBudgetOp[d]  = DInSquareOp[d] * ModBitsOp[d];
TwoTickPhaseSpaceBitBudgetOp[d_Integer] := TwoTickPhaseSpaceBitBudgetOp[d] = C0OverCSqOp[d] * SingleTickSpinorBitBudgetOp[d];

ShellSecondMomentTensorOp[d_Integer, shell_List] := Table[
  Sum[shell[[k, i]] * shell[[k, j]], {k, 1, Length[shell]}],
  {i, 1, d}, {j, 1, d}
];

ShellCrossFourthMomentOp[shell_List] :=
  Sum[shell[[k, 1]]^2 * shell[[k, 2]]^2, {k, 1, Length[shell]}];

ShellHasPlanarZeroCoordPred[shell_List] :=
  AnyTrue[shell, MemberQ[#, 0] &];
