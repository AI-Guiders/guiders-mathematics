(* ::Package:: *)
(* T <-> M Layer above Core CA: stream->lattice encoding, sensor readout, branchless 3D graph & Graphics3D *)

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

(* Encodes an arbitrary integer stream into a 3D VerletCell72 lattice of edge n *)
EncodeStreamToCarrier3DOp[d_Integer, stream_List, n_Integer] := Module[
  {len = Length[stream], nphi = NPhiOp[d], dphi = DeltaPhiDiscOp[d],
   nRing = NRingOp[d], kp = KPOp[d], muVac = MuVacOp[d], two = C0OverCSqOp[d]},
  Table[
    Module[{linIdx, v0, v1, mu0, mu1, p0, p1, p0Past, p1Past},
      linIdx = (iz - muVac) * n * n + (iy - muVac) * n + (ix - muVac);
      v0     = stream[[Mod[linIdx, len] + muVac]];
      v1     = stream[[Mod[linIdx + muVac, len] + muVac]];
      mu0    = EnforcePlanckCellFloorOp[d, Mod[(v0 + muVac) * nphi + linIdx, kp] + muVac];
      mu1    = EnforcePlanckCellFloorOp[d, Mod[(v1 + muVac) * nphi + two * linIdx, kp] + muVac];
      p0     = Mod[(v0 + muVac) * dphi + linIdx * SeamTicksOp[d], nRing];
      p1     = Mod[(v1 + muVac) * dphi + (linIdx + muVac) * SeamTicksOp[d], nRing];
      p0Past = Mod[p0 - dphi, nRing];
      p1Past = Mod[p1 - dphi, nRing];
      CreateVerletCell72Op[d, mu0, p0, p0Past, mu1, p1, p1Past]
    ],
    {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
  ]
];

(* Runs the pure 3D CA transition function StepCarrierForward3DOp for nSteps ticks *)
RunCarrier3DLabStepsOp[d_Integer, initGrid_List, nSteps_Integer] :=
  Nest[StepCarrierForward3DOp[d, #] &, initGrid, nSteps];

(* Pure branchless lattice readout without Select/filter:
   Every cell exists on the lattice; excitation indicator act(u) = IsPosNatOp[ne(u)] in {0, 1}.
   Node record layout: {1:Coord, 2:Mu0, 3:Phi0, 4:Mu1, 5:Phi1, 6:Gate0, 7:Ne, 8:Act, 9:Channel6, 10:RGB} *)
ExtractActiveLatticeGraph3DOp[d_Integer, grid_List] := Module[
  {n = Length[grid], two = C0OverCSqOp[d], nSq = NSquaresOp[d], dphi = DeltaPhiDiscOp[d],
   qSq = QSquareOp[d], muVac = MuVacOp[d],
   allCells, nTotal, fccPairs, nV, nE = 0, totalWeight = 0, nK3 = 0, nK4 = 0,
   nbUpper, channelCounts},
  allCells = Flatten[
    Table[
      Module[{cell = grid[[iz, iy, ix]], c3, s3, ne, act, ch, rgb},
        c3  = ProjCOp[d, cell[[2]]];
        s3  = ProjSOp[d, cell[[2]]];
        ne  = EnergyQuantumOp[d, cell[[7]]];
        act = IsPosNatOp[ne];
        ch  = Mod[Quotient[cell[[2]], dphi], nSq] + muVac;
        rgb = {Max[0, c3], Abs[s3], Max[0, -c3]} / qSq;
        {{iz, iy, ix}, cell[[1]], cell[[2]], cell[[4]], cell[[5]], cell[[7]], ne, act, ch, rgb}
      ],
      {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
    ],
    2
  ];
  nTotal = Length[allCells];
  nV     = Sum[allCells[[u, 8]], {u, 1, nTotal}];
  nbUpper = Table[
    Module[{cu = allCells[[u, 1]], list = {}},
      Do[
        Module[{du = cu - allCells[[v, 1]], isEdge},
          isEdge = EqIntIndicatorOp[du . du, two];
          list   = Join[list, ConstantArray[v, isEdge]]
        ],
        {v, u + 1, nTotal}
      ];
      list
    ],
    {u, 1, nTotal}
  ];
  fccPairs = Flatten[
    Table[
      Table[{u, v}, {v, nbUpper[[u]]}],
      {u, 1, nTotal}
    ],
    1
  ];
  Do[
    Module[{au = allCells[[u, 8]]},
      Do[
        Module[{av = allCells[[v, 8]], eAct},
          eAct = au * av;
          nE += eAct;
          totalWeight += eAct * (allCells[[u, 7]] + allCells[[v, 7]]);
          Do[
            Module[{aw = allCells[[w, 8]], tAct},
              tAct = eAct * aw;
              nK3 += tAct;
              Do[
                nK4 += tAct * allCells[[q, 8]],
                {q, Intersection[nbUpper[[u]], nbUpper[[v]], nbUpper[[w]]]}
              ]
            ],
            {w, Intersection[nbUpper[[u]], nbUpper[[v]]]}
          ]
        ],
        {v, nbUpper[[u]]}
      ]
    ],
    {u, 1, nTotal}
  ];
  channelCounts = Table[
    Sum[allCells[[u, 8]] * EqIntIndicatorOp[allCells[[u, 9]], c], {u, 1, nTotal}],
    {c, 1, nSq}
  ];
  <|"LatticeSize" -> nTotal,
    "NodeCount" -> nV,
    "EdgeCount" -> nE,
    "TriangleK3Count" -> nK3,
    "TetrahedraK4Count" -> nK4,
    "EulerChi2D" -> (nV - nE + nK3),
    "EulerChi3D" -> (nV - nE + nK3 - nK4),
    "MeanDegree" -> (two * nE) / Max[1, nV],
    "TotalEdgeWeight" -> totalWeight,
    "Channel6Spectrum" -> channelCounts,
    "Cells" -> allCells,
    "FccEdges" -> fccPairs|>
];

(* Pure branchless Graphics3D scene: every cell u has radius 0.08 * muVac + 0.05 * ne(u),
   and every FCC link (u, v) is drawn with opacity proportional to act(u) * act(v) *)
RenderLatticeGraphics3DOp[d_Integer, graphAssoc_Association] := Module[
  {cells = graphAssoc["Cells"], fccEdges = graphAssoc["FccEdges"], spherePrims, edgePrims},
  spherePrims = Table[
    {RGBColor @@ cells[[u, 10]],
     Sphere[cells[[u, 1]], 0.06 * cells[[u, 8]] * (MuVacOp[d] + cells[[u, 7]])]},
    {u, 1, Length[cells]}
  ];
  edgePrims = Table[
    Module[{u = fccEdges[[e, 1]], v = fccEdges[[e, 2]], eAct},
      eAct = cells[[u, 8]] * cells[[v, 8]];
      {Opacity[eAct], Line[{cells[[u, 1]], cells[[v, 1]]}]}
    ],
    {e, 1, Length[fccEdges]}
  ];
  Graphics3D[{spherePrims, edgePrims}, Boxed -> True]
];
