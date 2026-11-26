(* ::Package:: *)
(* Lab Layer above Core CA: Cuboctahedron coordination shells S_k(d), genetic code, 6 bio-metals *)

CuboctShellSizeOp[d_Integer, k_Integer] := Module[
  {isZero = IsZeroNatOp[k], shellK},
  shellK = (NNeighborsOp[d] - C0OverCSqOp[d]) * (k^C0OverCSqOp[d]) + C0OverCSqOp[d];
  BlendByIndicatorOp[isZero, MuVacOp[d], shellK]
];

CuboctBall2Op[d_Integer] :=
  MuVacOp[d] + CuboctShellSizeOp[d, MuVacOp[d]] + CuboctShellSizeOp[d, C0OverCSqOp[d]];

CuboctBallOp[d_Integer] :=
  CuboctBall2Op[d] + CuboctShellSizeOp[d, d];

HBondDivisorOp[d_Integer]      := KPOp[d] - MuVacOp[d];
CodonSpaceOp[d_Integer]        := DInSquareOp[d]^d;
AminoAcidClassesOp[d_Integer]  := NTrianglesOp[d] + NPhiOp[d];

WatsonCrickCompOp[d_Integer, b_Integer] := Mod[b + C0OverCSqOp[d], DInSquareOp[d]];
RumerStemInvolutionOp[d_Integer, s_Integer] :=
  (DInSquareOp[d]^C0OverCSqOp[d] - MuVacOp[d]) - s;

PlanetaryEnvWindowOp[d_Integer, mEnv_Integer] := NSquaresOp[d] + mEnv;
TransmembraneThicknessOp[d_Integer, m_Integer] := PlanetaryEnvWindowOp[d, m] * NNeighborsOp[d];

EnvSolventMatchGateOp[d_Integer, mGenome_Integer, mEnv_Integer] :=
  EqIntIndicatorOp[TransmembraneThicknessOp[d, mGenome], TransmembraneThicknessOp[d, mEnv]];

EnvMetabolicGateOp[d_Integer, phiEnv_Integer] :=
  IsPosNatOp[EnergyQuantumOp[d, HeisenbergDeadbandOp[d, phiEnv]]];

EnvFreeSurfaceGateOp[d_Integer, zInt_Integer] :=
  IsPosNatOp[NNeighborsOp[d] - zInt];

TissueEuler2DOp[complex_List] := complex[[1]] - complex[[2]] + complex[[3]];

ZygoteDormantComplexOp[d_Integer] := {MuVacOp[d], 0, 0, 0};
VesicleComplexOp[d_Integer] :=
  {NNeighborsOp[d], C0OverCSqOp[d] * NNeighborsOp[d], NTrianglesOp[d] + NSquaresOp[d], 0};

TubeComplexOp[d_Integer, r_Integer] :=
  {NSquaresOp[d] * r,
   C0OverCSqOp[d] * NTrianglesOp[d] * r - (NNeighborsOp[d] - C0OverCSqOp[d]),
   (NNeighborsOp[d] - C0OverCSqOp[d]) * (r - MuVacOp[d]),
   0};

CorticalColumnComplexOp[d_Integer, layers_Integer] :=
  {(NSquaresOp[d] + MuVacOp[d]) * layers,
   NNeighborsOp[d] * layers + d * DOutEdgeOp[d] * (layers - MuVacOp[d]),
   NSquaresOp[d] * layers + (NNeighborsOp[d] + C0OverCSqOp[d] * DOutEdgeOp[d]) * (layers - MuVacOp[d]),
   NSquaresOp[d] * (layers - MuVacOp[d])};

OntogenyGenesisInEnvOp[d_Integer, mGenome_Integer, sigmaAxis_Integer, sigmaIg_Integer,
  mEnv_Integer, phiEnv_Integer, rRings_Integer] :=
  BlendByIndicatorOp[
    EnvSolventMatchGateOp[d, mGenome, mEnv] * EnvMetabolicGateOp[d, phiEnv],
    BlendByIndicatorOp[
      sigmaAxis,
      BlendByIndicatorOp[
        sigmaIg,
        CorticalColumnComplexOp[d, NSquaresOp[d]],
        TubeComplexOp[d, rRings]
      ],
      VesicleComplexOp[d]
    ],
    ZygoteDormantComplexOp[d]
  ];

NuclearMagicLowOp[d_Integer, m_Integer] :=
  Quotient[m * (m + MuVacOp[d]) * (m + C0OverCSqOp[d]), d];

NuclearMagicHighOp[d_Integer, k_Integer] :=
  Quotient[k * (k^C0OverCSqOp[d] + DOutEdgeOp[d]), d];

ZHydrogenOp[d_Integer]  := MuVacOp[d];
ZHeliumOp[d_Integer]    := NuclearMagicLowOp[d, MuVacOp[d]];
ZLithiumOp[d_Integer]   := d;
ZBoronOp[d_Integer]     := DOutEdgeOp[d];
ZCarbonOp[d_Integer]    := NSquaresOp[d];
ZNitrogenOp[d_Integer]  := NTrianglesOp[d] - MuVacOp[d];
ZOxygenOp[d_Integer]    := NuclearMagicLowOp[d, C0OverCSqOp[d]];
ZNeonOp[d_Integer]      := NNeighborsOp[d] - C0OverCSqOp[d];
ZMagnesiumOp[d_Integer] := NNeighborsOp[d];
ZAluminumOp[d_Integer]  := NPhiOp[d];
ZSiliconOp[d_Integer]   := NPhiOp[d] + MuVacOp[d];
ZSulfurOp[d_Integer]    := C0OverCSqOp[d] * NTrianglesOp[d];
ZCalciumOp[d_Integer]   := NuclearMagicLowOp[d, d];
ZIronOp[d_Integer]      := NuclearMagicHighOp[d, DInSquareOp[d]] - C0OverCSqOp[d];
ZNickelOp[d_Integer]    := NuclearMagicHighOp[d, DInSquareOp[d]];
ZUraniumOp[d_Integer]   := CuboctShellSizeOp[d, d];
ZLeadTerminalOp[d_Integer] := NuclearMagicHighOp[d, NSquaresOp[d]];
NLeadTerminalOp[d_Integer] := NuclearMagicHighOp[d, NTrianglesOp[d] - MuVacOp[d]];

ClusterChargeSumOp[pairs_List] := Total[(#[[1]] * #[[2]]) & /@ pairs];

TetrapyrroleCageChargeOp[d_Integer]  := ClusterChargeSumOp[{{DInSquareOp[d], ZNitrogenOp[d]}}];
HexahistidineCageChargeOp[d_Integer] := ClusterChargeSumOp[{{NSquaresOp[d], ZNitrogenOp[d]}}];

DiscoverBioMetalInCageOp[maxZ_Integer, cageCharge_Integer, ligandCharge_Integer, targetShell_Integer] :=
  Select[Range[0, maxZ - 1], (# + cageCharge + ligandCharge == targetShell) &];

UniversalPhotonStepOp[d_Integer, mEnv_Integer] := Module[
  {mu = MuVacOp[d], hb = HBondDivisorOp[d]},
  EqIntIndicatorOp[mEnv, -mu] * ((NTrianglesOp[d] - mu) * hb) +
  EqIntIndicatorOp[mEnv, 0]   * (C0OverCSqOp[d] * hb) +
  EqIntIndicatorOp[mEnv, mu]  * hb
];
