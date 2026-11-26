(* ::Package:: *)
(* 3D Cuboctahedron Cellular Automaton: VerletCell72, 12-neighbor gate, IV(d) *)

EnergyQuantumOp[d_Integer, gateRaw_Integer] :=
  Quotient[Abs[Wrap512Op[d, gateRaw]], DeltaPhiDiscOp[d]];

HeisenbergDeadbandOp[d_Integer, gateRaw_Integer] := Module[
  {s = Wrap512Op[d, gateRaw], pass},
  pass = IsPosNatOp[EnergyQuantumOp[d, s]];
  pass * s
];

EnforcePlanckCellFloorOp[d_Integer, mu_Integer] :=
  Max[0, mu] + MuVacOp[d] * IsZeroNatOp[Max[0, mu]];

(* Pauli exclusion: fires Phi_Pauli(d) when both spinor lanes exceed K_P(d)^2 / 2 with |dphi| < DPhiDisc(d) *)
PauliExclusionGateOp[d_Integer, mu0_Integer, phi0_Integer, mu1_Integer, phi1_Integer, enabled_Integer] := Module[
  {rhoMinQ = Quotient[KPOp[d]^C0OverCSqOp[d], C0OverCSqOp[d]], d0, d1, par, fire},
  d0   = IsPosNatOp[mu0^C0OverCSqOp[d] - rhoMinQ];
  d1   = IsPosNatOp[mu1^C0OverCSqOp[d] - rhoMinQ];
  par  = IsZeroNatOp[EnergyQuantumOp[d, phi1 - phi0]];
  fire = IsPosNatOp[enabled] * d0 * d1 * par;
  PhiPauliOp[d] * fire
];

GateDenomOp[d_Integer, muC_Integer, zr_Integer] := Module[
  {r2 = (muC^C0OverCSqOp[d]) * QSquareOp[d], eps = EpsDenomOp[d], rawDen, isPos},
  rawDen = zr + r2 + eps;
  isPos  = IsPosNatOp[rawDen];
  BlendByIndicatorOp[isPos, rawDen, eps]
];

RawGatePhaseOp[d_Integer, muC_Integer, sumMu_Integer, zr_Integer, zi_Integer] := Module[
  {deficit, ziNonNeg, sgnZi, zDrive, safeDen, two = C0OverCSqOp[d], mu = MuVacOp[d]},
  deficit  = (sumMu * muC) * QSquareOp[d] - zr;
  ziNonNeg = mu - IsPosNatOp[-zi];
  sgnZi    = two * ziNonNeg - mu;
  zDrive   = zi + sgnZi * deficit;
  safeDen  = GateDenomOp[d, muC, zr];
  Wrap512Op[d, Floor[(KPOp[d] * zDrive) / safeDen]]
];

GatePhaseOp[d_Integer, muC_Integer, sumMu_Integer, zr_Integer, zi_Integer] :=
  HeisenbergDeadbandOp[d, RawGatePhaseOp[d, muC, sumMu, zr, zi]];

Compute12NeighborGateFromCellsOp[d_Integer, muC_Integer, phiC_Integer, neighbors12_List] := Module[
  {sumMu = 0, zrRaw = 0, ziRaw = 0, mN, pN, dphi},
  Do[
    mN = neighbors12[[k, 1]];
    pN = neighbors12[[k, 2]];
    dphi = pN - phiC;
    sumMu += mN;
    zrRaw += mN * ProjCOp[d, dphi];
    ziRaw += mN * ProjSOp[d, dphi],
    {k, 1, Length[neighbors12]}
  ];
  GatePhaseOp[d, muC, sumMu, muC * zrRaw, muC * ziRaw]
];

SelectLeadingGateOp[d_Integer, g0_Integer, g1_Integer] := Module[
  {a0 = Abs[Wrap512Op[d, g0]], a1 = Abs[Wrap512Op[d, g1]], take1},
  take1 = IsPosNatOp[a1 - a0];
  BlendByIndicatorOp[take1, g1, g0]
];

ComputeTwoLaneCellGateOp[d_Integer,
  mu0_Integer, phi0_Integer, mu1_Integer, phi1_Integer,
  nbLane0_List, nbLane1_List, pauliEnabled_Integer] := Module[
  {g0, g1, lead, pauli},
  g0    = Compute12NeighborGateFromCellsOp[d, mu0, phi0, nbLane0];
  g1    = Compute12NeighborGateFromCellsOp[d, mu1, phi1, nbLane1];
  lead  = SelectLeadingGateOp[d, g0, g1];
  pauli = PauliExclusionGateOp[d, mu0, phi0, mu1, phi1, pauliEnabled];
  Mod[lead + pauli, NRingOp[d]]
];

(* VerletCell72 layout: {mu0, phi0Curr, phi0Past, mu1, phi1Curr, phi1Past, gate0, gate1}.
   Axiom A5 (gapless occupancy): every cell in the lattice is filled (mu0 >= mu_vac(d), mu1 >= mu_vac(d)). *)
CreateVerletCell72Op[d_Integer, mu0_Integer, p0Curr_Integer, p0Past_Integer, mu1_Integer, p1Curr_Integer, p1Past_Integer] :=
  {EnforcePlanckCellFloorOp[d, Mod[mu0, NRingOp[d]]],
   Mod[p0Curr, NRingOp[d]],
   Mod[p0Past, NRingOp[d]],
   EnforcePlanckCellFloorOp[d, Mod[mu1, NRingOp[d]]],
   Mod[p1Curr, NRingOp[d]],
   Mod[p1Past, NRingOp[d]],
   0, 0};

StepCellForwardOp[d_Integer, cell_List, gate0_Integer, gate1_Integer] := Module[
  {mu0 = cell[[1]], p0Curr = cell[[2]], p0Past = cell[[3]],
   mu1 = cell[[4]], p1Curr = cell[[5]], p1Past = cell[[6]],
   nRing = NRingOp[d], two = C0OverCSqOp[d], p0Next, p1Next},
  p0Next = Mod[two * p0Curr - p0Past + gate0, nRing];
  p1Next = Mod[two * p1Curr - p1Past + gate1, nRing];
  {mu0, p0Next, p0Curr, mu1, p1Next, p1Curr, Mod[gate0, nRing], Mod[gate1, nRing]}
];

StepCellReverseOp[d_Integer, steppedCell_List, gate0AtT_Integer, gate1AtT_Integer] := Module[
  {mu0 = steppedCell[[1]], p0Next = steppedCell[[2]], p0Curr = steppedCell[[3]],
   mu1 = steppedCell[[4]], p1Next = steppedCell[[5]], p1Curr = steppedCell[[6]],
   nRing = NRingOp[d], two = C0OverCSqOp[d], p0PastRec, p1PastRec},
  p0PastRec = Mod[two * p0Curr - p0Next + gate0AtT, nRing];
  p1PastRec = Mod[two * p1Curr - p1Next + gate1AtT, nRing];
  {mu0, p0Curr, p0PastRec, mu1, p1Curr, p1PastRec, Mod[gate0AtT, nRing], Mod[gate1AtT, nRing]}
];

(* Branchless bitwise integer square root floor(sqrt(q)) on [0, K_P(d)^2] *)
IsqrtOp[d_Integer, q_Integer] := Module[{r = 0, bit, cand, fits},
  Do[
    bit  = C0OverCSqOp[d] ^ k;
    cand = r + bit;
    fits = IsZeroNatOp[Quotient[cand^C0OverCSqOp[d], Max[0, q] + MuVacOp[d]]];
    r    = r + fits * bit,
    {k, FracBitsOp[d], 0, -1}
  ];
  r
];

ComputeSingleDirectedEdgeFluxOp[d_Integer, muX_Integer, phiX_Integer, muY_Integer, phiY_Integer] := Module[
  {dphiXY = phiX - phiY, jRaw, qX, qY, qMax, halfLinks, capPos, capNeg, isNonNeg},
  jRaw = Floor[(muX * muY * ProjSOp[d, dphiXY]) / (NNeighborsOp[d] * QSquareOp[d])];
  qX   = muX^C0OverCSqOp[d];
  qY   = muY^C0OverCSqOp[d];
  qMax = KPOp[d]^C0OverCSqOp[d];
  halfLinks = NSquaresOp[d];
  capPos = Min[Quotient[qX - MuVacOp[d], halfLinks], Quotient[qMax - qY, halfLinks]];
  capNeg = Min[Quotient[qY - MuVacOp[d], halfLinks], Quotient[qMax - qX, halfLinks]];
  isNonNeg = MuVacOp[d] - IsPosNatOp[-jRaw];
  BlendByIndicatorOp[isNonNeg, Min[Max[0, jRaw], capPos], -Min[Max[0, -jRaw], capNeg]]
];

(* Conservative 12-edge mass transport preserving gapless occupancy mu_next >= mu_vac(d) *)
Compute12EdgeMassFluxStepAtCellOp[d_Integer, muC_Integer, phiC_Integer, posNegPairs6_List] := Module[
  {divJ = 0, mPos, pPos, mNeg, pNeg, jOut, jIn, qCurr, qNext, muNext},
  Do[
    mPos = posNegPairs6[[k, 1]]; pPos = posNegPairs6[[k, 2]];
    mNeg = posNegPairs6[[k, 3]]; pNeg = posNegPairs6[[k, 4]];
    jOut = ComputeSingleDirectedEdgeFluxOp[d, muC, phiC, mPos, pPos];
    jIn  = ComputeSingleDirectedEdgeFluxOp[d, mNeg, pNeg, muC, phiC];
    divJ += (jOut - jIn),
    {k, 1, Length[posNegPairs6]}
  ];
  qCurr  = muC^C0OverCSqOp[d];
  qNext  = qCurr - divJ;
  muNext = EnforcePlanckCellFloorOp[d, IsqrtOp[d, qNext]];
  {muNext, divJ}
];

ReadCellSensorsOp[d_Integer, cell_List] := Module[
  {mu0 = cell[[1]], phi0 = cell[[2]], mu1 = cell[[4]], gate0 = cell[[7]],
   c3, s3, rCh, gCh, bCh, sound, rxn, ne, shellZ, centerNm, wavelengthNm},
  c3 = ProjCOp[d, phi0];
  s3 = ProjSOp[d, phi0];
  rCh = Max[0, c3];
  gCh = Abs[s3];
  bCh = Max[0, -c3];
  sound = mu0 * mu1 * (QSquareOp[d] - c3);
  rxn   = mu0 * mu1 * s3;
  ne    = EnergyQuantumOp[d, gate0];
  shellZ = mu0 + mu1;
  centerNm = NRingOp[d] + KPOp[d] + d;
  wavelengthNm = centerNm + Floor[(CollapsePhiOp[d] * c3) / QSquareOp[d]];
  <|"ColorR" -> rCh, "ColorG" -> gCh, "ColorB" -> bCh,
    "AcousticWave" -> sound, "ReactionCurrent" -> rxn,
    "EnergyQuantaNe" -> ne, "ShellZ" -> shellZ, "WavelengthNm" -> wavelengthNm|>
];

Iv263InBallIndOp[d_Integer, dz_Integer, dy_Integer, dx_Integer] := Module[
  {r2 = dz^C0OverCSqOp[d] + dy^C0OverCSqOp[d] + dx^C0OverCSqOp[d],
   maxR2 = NSquaresOp[d]^C0OverCSqOp[d]},
  IsZeroNatOp[Quotient[r2, maxR2 + MuVacOp[d]]]
];

Iv263MuAtOp[d_Integer, dz_Integer, dy_Integer, dx_Integer] :=
  MuVacOp[d] + NNeighborsOp[d] * Iv263InBallIndOp[d, dz, dy, dx];

Iv263PhiAtOp[d_Integer, dz_Integer, dy_Integer, dx_Integer, lane_Integer, centerX_Integer: 0] := Module[
  {mu = MuVacOp[d], two = C0OverCSqOp[d], four = DInSquareOp[d],
   cx, zNonNeg, yNonNeg, xNonNeg, isLane1, yEff, xEff, shift, mask8, bit, bg, shiftTicks, inBall},
  cx = BlendByIndicatorOp[IsZeroNatOp[centerX], NNeighborsOp[d], centerX];
  zNonNeg = mu - IsPosNatOp[-dz];
  yNonNeg = mu - IsPosNatOp[-dy];
  xNonNeg = mu - IsPosNatOp[-dx];
  isLane1 = Mod[lane, two];
  yEff = BlendByIndicatorOp[isLane1, mu - yNonNeg, yNonNeg];
  xEff = BlendByIndicatorOp[isLane1, mu - xNonNeg, xNonNeg];
  shift = four * zNonNeg + two * yEff + mu * xEff;
  mask8 = Mod[IvSeedOp[d], PhiPauliOp[d]];
  bit = Mod[Quotient[mask8, two^shift], two];
  bg = Mod[Mod[dx + cx, NPhiOp[d]] * DeltaPhiDiscOp[d], NRingOp[d]];
  shiftTicks = d * DeltaPhiDiscOp[d];
  inBall = Iv263InBallIndOp[d, dz, dy, dx];
  Mod[bg + inBall * bit * shiftTicks, NRingOp[d]]
];

Compute12NeighborHolonomyGateAtIvOp[d_Integer, dz_Integer, dy_Integer, dx_Integer, lane_Integer] := Module[
  {muC, phiC, offsets = O2CuboctahedronVerticesOp[d], nbs},
  muC  = Iv263MuAtOp[d, dz, dy, dx];
  phiC = Iv263PhiAtOp[d, dz, dy, dx, lane];
  nbs  = Table[
    {Iv263MuAtOp[d, dz + offsets[[k, 1]], dy + offsets[[k, 2]], dx + offsets[[k, 3]]],
     Iv263PhiAtOp[d, dz + offsets[[k, 1]], dy + offsets[[k, 2]], dx + offsets[[k, 3]], lane]},
    {k, 1, Length[offsets]}
  ];
  Compute12NeighborGateFromCellsOp[d, muC, phiC, nbs]
];

(* Full 3D lattice initialized from seed IV(d) = 263: every cell is filled (mu = 13 in ball, mu = 1 in vacuum) *)
CreateIvCarrier3DOp[d_Integer, n_Integer] := Module[
  {c = Quotient[n, C0OverCSqOp[d]], one = MuVacOp[d]},
  Table[
    Module[{dz = iz - one - c, dy = iy - one - c, dx = ix - one - c, mu, p0, p1},
      mu = Iv263MuAtOp[d, dz, dy, dx];
      p0 = Iv263PhiAtOp[d, dz, dy, dx, 0, c];
      p1 = Iv263PhiAtOp[d, dz, dy, dx, one, c];
      CreateVerletCell72Op[d, mu, p0, p0, mu, p1, p1]
    ],
    {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
  ]
];

AllCellsFilledPred[d_Integer, grid_List] := Module[{muMin = MuVacOp[d]},
  AllTrue[Flatten[grid, 2], (#[[1]] >= muMin && #[[4]] >= muMin) &]
];

CreateReactionCarrier3DOp[d_Integer, n_Integer, phiEnv_Integer] := Module[
  {c = Quotient[n, C0OverCSqOp[d]], mu = NPhiOp[d], two = C0OverCSqOp[d], four = DInSquareOp[d], one = MuVacOp[d]},
  Table[
    Module[{dz = iz - one - c, dy = iy - one - c, dx = ix - one - c, isCenter, lin, isPosHalf, envPhase, p},
      isCenter = IsZeroNatOp[Abs[dz] + Abs[dy] + Abs[dx]];
      lin = four * dz + two * dy + dx;
      isPosHalf = one - IsPosNatOp[-lin];
      envPhase = BlendByIndicatorOp[isPosHalf, Mod[phiEnv, NRingOp[d]], Mod[-phiEnv, NRingOp[d]]];
      p = BlendByIndicatorOp[isCenter, 0, envPhase];
      CreateVerletCell72Op[d, mu, p, p, mu, p, p]
    ],
    {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
  ]
];

StepCarrierGrid3DOp[d_Integer, grid_List, usePast_Integer, isReverse_Integer] := Module[
  {n = Length[grid], offsets = O2CuboctahedronVerticesOp[d], pIdx0, pIdx1},
  pIdx0 = BlendByIndicatorOp[usePast, 3, 2];
  pIdx1 = BlendByIndicatorOp[usePast, 6, 5];
  Table[
    Module[{cell = grid[[iz, iy, ix]], nb0, nb1, gate, fwdCell, revCell},
      nb0 = Table[
        Module[{nz = Mod[iz - 1 + offsets[[k, 1]], n] + 1,
                ny = Mod[iy - 1 + offsets[[k, 2]], n] + 1,
                nx = Mod[ix - 1 + offsets[[k, 3]], n] + 1, nbCell},
          nbCell = grid[[nz, ny, nx]];
          {nbCell[[1]], nbCell[[pIdx0]]}
        ],
        {k, 1, Length[offsets]}
      ];
      nb1 = Table[
        Module[{nz = Mod[iz - 1 + offsets[[k, 1]], n] + 1,
                ny = Mod[iy - 1 + offsets[[k, 2]], n] + 1,
                nx = Mod[ix - 1 + offsets[[k, 3]], n] + 1, nbCell},
          nbCell = grid[[nz, ny, nx]];
          {nbCell[[4]], nbCell[[pIdx1]]}
        ],
        {k, 1, Length[offsets]}
      ];
      gate = ComputeTwoLaneCellGateOp[d, cell[[1]], cell[[pIdx0]], cell[[4]], cell[[pIdx1]], nb0, nb1, MuVacOp[d]];
      fwdCell = StepCellForwardOp[d, cell, gate, gate];
      revCell = StepCellReverseOp[d, cell, gate, gate];
      BlendByIndicatorOp[isReverse, revCell, fwdCell]
    ],
    {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
  ]
];

StepCarrierForward3DOp[d_Integer, grid_List]  := StepCarrierGrid3DOp[d, grid, 0, 0];
StepCarrierBackward3DOp[d_Integer, grid_List] := StepCarrierGrid3DOp[d, grid, MuVacOp[d], MuVacOp[d]];

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

  (EnforcePlanckCellFloorOp[d, 0] == MuVacOp[d]) &&
  AllCellsFilledPred[d, ivGrid] &&
  AllCellsFilledPred[d, bzGrid1] &&
  (g0 == 0) &&
  (gWallPos + gWallNeg == -MuVacOp[d]) &&
  (inert1[[cIdx, cIdx, cIdx, 2]] == 0) &&
  (bzGrid1[[cIdx, cIdx, cIdx, 2]] == NRingOp[d] - DOutEdgeOp[d] * Quotient[QSquareOp[d], DInSquareOp[d]]) &&
  (bzRev[[All, All, All, 1 ;; 6]] == bzGrid0[[All, All, All, 1 ;; 6]])
];
