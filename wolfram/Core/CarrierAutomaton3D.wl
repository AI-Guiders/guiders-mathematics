(* ::Package:: *)
(* Pure M-Layer Core: 3D Cuboctahedron Cellular Automaton (VerletCell72, 12-neighbor gate, 12-edge mass flux) *)

EnergyQuantumOp[d_Integer, gateRaw_Integer] :=
  Quotient[Abs[Wrap512Op[d, gateRaw]], DeltaPhiDiscOp[d]];

HeisenbergDeadbandOp[d_Integer, gateRaw_Integer] := Module[
  {s = Wrap512Op[d, gateRaw], pass},
  pass = IsPosNatOp[EnergyQuantumOp[d, s]];
  pass * s
];

EnforcePlanckCellFloorOp[d_Integer, mu_Integer] :=
  Max[0, mu] + MuVacOp[d] * IsZeroNatOp[Max[0, mu]];

RingVisibilityL1NormOp[d_Integer, mu_Integer, phi_Integer] :=
  mu * (Abs[ProjCOp[d, phi]] + Abs[ProjSOp[d, phi]]);

PauliExclusionGateOp[d_Integer, mu0_Integer, phi0_Integer, mu1_Integer, phi1_Integer, enabled_Integer] := Module[
  {two = C0OverCSqOp[d], kp = KPOp[d], rhoMinQ, d0, d1, par, fire},
  rhoMinQ = Quotient[kp^two, two];
  d0   = IsPosNatOp[mu0^two - rhoMinQ];
  d1   = IsPosNatOp[mu1^two - rhoMinQ];
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
  {lutC = RingFaceLutCOp[d], lutS = RingFaceLutSOp[d], nRing = NRingOp[d], mu = MuVacOp[d],
   sumMu = 0, zrRaw = 0, ziRaw = 0, mN, idx},
  Do[
    mN = neighbors12[[k, 1]];
    idx = Mod[neighbors12[[k, 2]] - phiC, nRing] + mu;
    sumMu += mN;
    zrRaw += mN * lutC[[idx]];
    ziRaw += mN * lutS[[idx]],
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

CreateVerletCell72Op[d_Integer, mu0_Integer, p0Curr_Integer, p0Past_Integer, mu1_Integer, p1Curr_Integer, p1Past_Integer] := Module[
  {nRing = NRingOp[d]},
  {EnforcePlanckCellFloorOp[d, Mod[mu0, nRing]],
   Mod[p0Curr, nRing],
   Mod[p0Past, nRing],
   EnforcePlanckCellFloorOp[d, Mod[mu1, nRing]],
   Mod[p1Curr, nRing],
   Mod[p1Past, nRing],
   0, 0}
];

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

IsqrtLutOp[d_Integer] := IsqrtLutOp[d] = Module[
  {qMax = KPOp[d]^C0OverCSqOp[d]},
  Table[Floor[Sqrt[q]], {q, 0, qMax}]
];

IsqrtOp[d_Integer, q_Integer] := Module[
  {qClamped = Min[Max[0, q], KPOp[d]^C0OverCSqOp[d]]},
  IsqrtLutOp[d][[qClamped + MuVacOp[d]]]
];

ComputeSingleDirectedEdgeFluxOp[d_Integer, muX_Integer, phiX_Integer, muY_Integer, phiY_Integer] := Module[
  {lutS = RingFaceLutSOp[d], nRing = NRingOp[d], muVac = MuVacOp[d], two = C0OverCSqOp[d],
   denom = NNeighborsOp[d] * QSquareOp[d], qMax = KPOp[d]^C0OverCSqOp[d], halfLinks = NSquaresOp[d],
   sVal, jRaw, qX, qY, capPos, capNeg, isNonNeg},
  sVal = lutS[[Mod[phiX - phiY, nRing] + muVac]];
  jRaw = Floor[(muX * muY * sVal) / denom];
  qX   = muX^two;
  qY   = muY^two;
  capPos = Min[Quotient[qX - muVac, halfLinks], Quotient[qMax - qY, halfLinks]];
  capNeg = Min[Quotient[qY - muVac, halfLinks], Quotient[qMax - qX, halfLinks]];
  isNonNeg = muVac - IsPosNatOp[-jRaw];
  BlendByIndicatorOp[isNonNeg, Min[Max[0, jRaw], capPos], -Min[Max[0, -jRaw], capNeg]]
];

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
  {n = Length[grid], offsets = O2CuboctahedronVerticesOp[d], nOff, muVac = MuVacOp[d], pIdx0, pIdx1},
  nOff  = Length[offsets];
  pIdx0 = BlendByIndicatorOp[usePast, 3, 2];
  pIdx1 = BlendByIndicatorOp[usePast, 6, 5];
  Table[
    Module[{cell = grid[[iz, iy, ix]], nb0, nb1, gate, fwdCell, revCell},
      nb0 = Table[
        Module[{nbCell = grid[[Mod[iz - 1 + offsets[[k, 1]], n] + 1,
                               Mod[iy - 1 + offsets[[k, 2]], n] + 1,
                               Mod[ix - 1 + offsets[[k, 3]], n] + 1]]},
          {nbCell[[1]], nbCell[[pIdx0]]}
        ],
        {k, 1, nOff}
      ];
      nb1 = Table[
        Module[{nbCell = grid[[Mod[iz - 1 + offsets[[k, 1]], n] + 1,
                               Mod[iy - 1 + offsets[[k, 2]], n] + 1,
                               Mod[ix - 1 + offsets[[k, 3]], n] + 1]]},
          {nbCell[[4]], nbCell[[pIdx1]]}
        ],
        {k, 1, nOff}
      ];
      gate = ComputeTwoLaneCellGateOp[d, cell[[1]], cell[[pIdx0]], cell[[4]], cell[[pIdx1]], nb0, nb1, muVac];
      fwdCell = StepCellForwardOp[d, cell, gate, gate];
      revCell = StepCellReverseOp[d, cell, gate, gate];
      BlendByIndicatorOp[isReverse, revCell, fwdCell]
    ],
    {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
  ]
];

StepCarrierForward3DOp[d_Integer, grid_List]  := StepCarrierGrid3DOp[d, grid, 0, 0];
StepCarrierBackward3DOp[d_Integer, grid_List] := StepCarrierGrid3DOp[d, grid, MuVacOp[d], MuVacOp[d]];
