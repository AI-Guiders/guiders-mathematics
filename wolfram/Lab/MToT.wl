(* ::Package:: *)
(* T <-> M Layer above Core CA: stream->lattice encoding, sensor readout, 3D graph & Graphics3D extraction *)

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

(* Runs the pure 3D CA (phase Verlet + 12-edge mass flux) for nSteps ticks *)
RunCarrier3DLabStepsOp[d_Integer, initGrid_List, nSteps_Integer] :=
  Nest[StepCarrierWithMassFlux3DOp[d, #] &, initGrid, nSteps];

(* Node record layout: {1:Coord, 2:Mu0, 3:Phi0, 4:Mu1, 5:Phi1, 6:Gate0, 7:EnergyQuanta, 8:Channel6, 9:RGB} *)
ExtractActiveLatticeGraph3DOp[d_Integer, grid_List] := Module[
  {n = Length[grid], two = C0OverCSqOp[d], nSq = NSquaresOp[d], dphi = DeltaPhiDiscOp[d],
   qSq = QSquareOp[d], muVac = MuVacOp[d], nphi = NPhiOp[d],
   allCells, activeNodes, nV, nbList, edges = {}, edgeWeights = {}, tris = {}, tets = {},
   channelCounts},
  allCells = Flatten[
    Table[
      Module[{cell = grid[[iz, iy, ix]], c3, s3, ne, ch, rgb},
        c3  = ProjCOp[d, cell[[2]]];
        s3  = ProjSOp[d, cell[[2]]];
        ne  = EnergyQuantumOp[d, cell[[7]]];
        ch  = Mod[Quotient[cell[[2]], dphi], nSq] + muVac;
        rgb = {Max[0, c3], Abs[s3], Max[0, -c3]} / qSq;
        {{iz, iy, ix}, cell[[1]], cell[[2]], cell[[4]], cell[[5]], cell[[7]], ne, ch, rgb}
      ],
      {iz, 1, n}, {iy, 1, n}, {ix, 1, n}
    ],
    2
  ];
  activeNodes = Select[allCells, (#[[7]] > 0 || (#[[2]] + #[[4]]) >= two * nphi) &];
  nV = Length[activeNodes];
  nbList = Table[
    Select[Range[u + 1, nV],
      Module[{du = activeNodes[[u, 1]] - activeNodes[[#, 1]], flux, dGate},
        ((du . du) == two) &&
        (
          flux  = Abs[ComputeSingleDirectedEdgeFluxOp[d,
                    activeNodes[[u, 2]], activeNodes[[u, 3]],
                    activeNodes[[#, 2]], activeNodes[[#, 3]]]];
          dGate = EnergyQuantumOp[d, activeNodes[[u, 3]] - activeNodes[[#, 3]]];
          (flux + dGate) > 0
        )
      ] &
    ],
    {u, 1, nV}
  ];
  Do[
    Do[
      AppendTo[edges, {u, v}];
      AppendTo[edgeWeights,
        Abs[ComputeSingleDirectedEdgeFluxOp[d,
              activeNodes[[u, 2]], activeNodes[[u, 3]],
              activeNodes[[v, 2]], activeNodes[[v, 3]]]] +
        EnergyQuantumOp[d, activeNodes[[u, 3]] - activeNodes[[v, 3]]]
      ];
      Do[
        AppendTo[tris, {u, v, w}];
        Do[
          AppendTo[tets, {u, v, w, q}],
          {q, Intersection[nbList[[u]], nbList[[v]], nbList[[w]]]}
        ],
        {w, Intersection[nbList[[u]], nbList[[v]]]}
      ],
      {v, nbList[[u]]}
    ],
    {u, 1, nV}
  ];
  channelCounts = Table[Length[Select[activeNodes, (#[[8]] == c) &]], {c, 1, nSq}];
  <|"NodeCount" -> nV,
    "EdgeCount" -> Length[edges],
    "TriangleK3Count" -> Length[tris],
    "TetrahedraK4Count" -> Length[tets],
    "EulerChi2D" -> (nV - Length[edges] + Length[tris]),
    "MeanDegree" -> (two * Length[edges]) / Max[1, nV],
    "Channel6Spectrum" -> channelCounts,
    "Nodes" -> activeNodes,
    "Edges" -> edges,
    "EdgeWeights" -> edgeWeights,
    "TrianglesK3" -> tris,
    "TetrahedraK4" -> tets|>
];

(* Builds a symbolic Wolfram Language Graphics3D scene from the extracted 3D lattice graph *)
RenderLatticeGraphics3DOp[d_Integer, graphAssoc_Association] := Module[
  {nodes = graphAssoc["Nodes"], edges = graphAssoc["Edges"], spherePrims, edgePrims},
  spherePrims = Table[
    {RGBColor @@ nodes[[u, 9]],
     Sphere[nodes[[u, 1]], 0.28]},
    {u, 1, Length[nodes]}
  ];
  edgePrims = Table[
    Line[{nodes[[edges[[e, 1]], 1]], nodes[[edges[[e, 2]], 1]]}],
    {e, 1, Length[edges]}
  ];
  Graphics3D[{spherePrims, edgePrims}, Boxed -> True]
];
