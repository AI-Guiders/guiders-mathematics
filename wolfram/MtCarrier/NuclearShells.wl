(* ::Package:: *)
(* Generative Periodic Table (Periods 1..8), 400-junction ATG protection, DNA->tissue compiler *)

ChannelOrbitalCountOp[d_Integer, ell_Integer] :=
  C0OverCSqOp[d] * ell + MuVacOp[d];

SubshellCapacityOp[d_Integer, ell_Integer] :=
  C0OverCSqOp[d] * ChannelOrbitalCountOp[d, ell];

MadelungPeriodLengthOp[d_Integer, p_Integer] := Module[
  {two = C0OverCSqOp[d]},
  two * (Quotient[p + two, two] ^ two)
];

NobleGasZOp[d_Integer, p_Integer] :=
  Total[Table[MadelungPeriodLengthOp[d, q], {q, 1, p}]];

ValenceElectronsInPeriodOp[d_Integer, p_Integer, z_Integer] :=
  z - NobleGasZOp[d, p - MuVacOp[d]];

CentralGreenPotentialJumpOp[d_Integer] :=
  Quotient[d, (MuVacOp[d]^d - 0^d)];

ZOganessonOp[d_Integer]               := NobleGasZOp[d, NTrianglesOp[d] - MuVacOp[d]];
ZIslandFlOp[d_Integer]                := NLeadTerminalOp[d] - NNeighborsOp[d];
ZIslandUbnOp[d_Integer]               := NLeadTerminalOp[d] - NSquaresOp[d];
ZCoreDiracDivingOp[d_Integer]         := d * DeltaPhiDiscOp[d];
ZFissionCeilingOp[d_Integer]          := NLeadTerminalOp[d];
ZFeynmanPointLimitOp[d_Integer]       := AlphaFloorIntOp[d];
Period8RemainingElementsOp[d_Integer] := ZFeynmanPointLimitOp[d] - ZOganessonOp[d];
ZQedFiniteNucleusLowOp[d_Integer]     := ZFeynmanPointLimitOp[d] + DOutEdgeOp[d] * (NTrianglesOp[d] - MuVacOp[d]);
ZQedFiniteNucleusHighOp[d_Integer]    := ZFeynmanPointLimitOp[d] + NSquaresOp[d]^C0OverCSqOp[d];
NIslandOfStabilityOp[d_Integer]       := NuclearMagicHighOp[d, NTrianglesOp[d]];

CodonTripletIndexOp[d_Integer, b1_Integer, b2_Integer, b3_Integer] :=
  (DInSquareOp[d]^C0OverCSqOp[d]) * b1 + DInSquareOp[d] * b2 + b3;

CodonBase1Op[d_Integer, cod_Integer] :=
  Mod[Quotient[cod, DInSquareOp[d]^C0OverCSqOp[d]], DInSquareOp[d]];

CodonBase2Op[d_Integer, cod_Integer] :=
  Mod[Quotient[cod, DInSquareOp[d]], DInSquareOp[d]];

CodonBase3Op[d_Integer, cod_Integer] :=
  Mod[cod, DInSquareOp[d]];

StartCodonAtgOp[d_Integer] :=
  CodonTripletIndexOp[d, C0OverCSqOp[d], 0, d];

(* Canonical triplet codon assignment: b_3 % 2 == mu_vac(d), preventing out-of-frame ATG at shifts +1, +2 *)
AminoAcidToTripletCodonOp[d_Integer, a_Integer] := Module[
  {four = DInSquareOp[d], sixteen = DInSquareOp[d]^C0OverCSqOp[d]},
  four * Mod[a, sixteen] + (MuVacOp[d] + C0OverCSqOp[d] * Quotient[a, sixteen])
];

TripletCodonToAminoAcidOp[d_Integer, cod_Integer] := Module[
  {four = DInSquareOp[d], sixteen = DInSquareOp[d]^C0OverCSqOp[d]},
  Quotient[cod, four] + sixteen * Quotient[Mod[cod, four], d]
];

IsJunctionFrameProtectedPred[d_Integer, cod1_Integer, cod2_Integer] := Module[
  {atg = StartCodonAtgOp[d], shift1, shift2},
  shift1 = CodonTripletIndexOp[d, CodonBase2Op[d, cod1], CodonBase3Op[d, cod1], CodonBase1Op[d, cod2]];
  shift2 = CodonTripletIndexOp[d, CodonBase3Op[d, cod1], CodonBase1Op[d, cod2], CodonBase2Op[d, cod2]];
  (shift1 != atg) && (shift2 != atg)
];

InverseTargetSigmaAxisOp[d_Integer, chiTarget_Integer] :=
  MuVacOp[d] - EqIntIndicatorOp[chiTarget, C0OverCSqOp[d]];

InverseTargetSigmaIgOp[d_Integer, chiTarget_Integer] :=
  EqIntIndicatorOp[chiTarget, DeltaPhiDiscOp[d]];

GrooveSynthOp[d_Integer, b_Integer] := Module[
  {mu = MuVacOp[d], two = C0OverCSqOp[d]},
  EqIntIndicatorOp[b, 0]   * d +
  EqIntIndicatorOp[b, mu]  * (AminoAcidClassesOp[d] - two) +
  EqIntIndicatorOp[b, two] * (d * DOutEdgeOp[d]) +
  EqIntIndicatorOp[b, d]   * two
];

GrooveReadOp[d_Integer, a_Integer] := Module[
  {mu = MuVacOp[d], two = C0OverCSqOp[d]},
  EqIntIndicatorOp[a, d] * 0 +
  EqIntIndicatorOp[a, AminoAcidClassesOp[d] - two] * mu +
  EqIntIndicatorOp[a, d * DOutEdgeOp[d]] * two +
  EqIntIndicatorOp[a, two] * d
];

CompilePhenotypeToDnaSequenceOp[d_Integer, chiTarget_Integer] := Module[
  {codAxis, codIg},
  codAxis = AminoAcidToTripletCodonOp[d, GrooveSynthOp[d, InverseTargetSigmaAxisOp[d, chiTarget]]];
  codIg   = AminoAcidToTripletCodonOp[d, GrooveSynthOp[d, InverseTargetSigmaIgOp[d, chiTarget]]];
  {CodonBase1Op[d, codAxis], CodonBase2Op[d, codAxis], CodonBase3Op[d, codAxis],
   CodonBase1Op[d, codIg],   CodonBase2Op[d, codIg],   CodonBase3Op[d, codIg]}
];

TranslateTwoCodonDnaToProteinOp[d_Integer, dna_List] :=
  {TripletCodonToAminoAcidOp[d, CodonTripletIndexOp[d, dna[[1]], dna[[2]], dna[[3]]]],
   TripletCodonToAminoAcidOp[d, CodonTripletIndexOp[d, dna[[4]], dna[[5]], dna[[6]]]]};

GrowTissueFromDnaSequenceInEnvOp[d_Integer, dna_List,
  mGenome_Integer, mEnv_Integer, phiEnv_Integer, rRings_Integer] := Module[
  {prot = TranslateTwoCodonDnaToProteinOp[d, dna]},
  OntogenyGenesisInEnvOp[
    d, mGenome,
    GrooveReadOp[d, prot[[1]]],
    GrooveReadOp[d, prot[[2]]],
    mEnv, phiEnv, rRings
  ]
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
  (TissueEuler2DOp[GrowTissueFromDnaSequenceInEnvOp[d, CompilePhenotypeToDnaSequenceOp[d, dphi], 0, 0, dphi, d]] == dphi)
];
