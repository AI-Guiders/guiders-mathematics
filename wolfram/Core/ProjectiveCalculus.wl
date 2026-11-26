(* ::Package:: *)
(* Pure Mathematical Projective Manifold & CAS Calculus over Z_{2^m} for arbitrary rank m -> Infinity,
   arbitrary spatial dimension d >= 2, and arbitrary integration dimension nDim >= 1 (zero physical locks) *)

ProjectiveRingSizeOp[m_Integer]     := 2 ^ m;
ProjectiveHalfTurnOp[m_Integer]     := 2 ^ (m - 1);
ProjectiveQuadrantSizeOp[m_Integer] := 2 ^ (m - 2);

ProjectiveWrapOp[m_Integer, phi_Integer] := Module[
  {nRing = ProjectiveRingSizeOp[m], half = ProjectiveHalfTurnOp[m], u, q},
  u = Mod[phi, nRing];
  q = Quotient[u + (half - 1), nRing];
  u - nRing * q
];

ProjectiveCharCOp[m_Integer, phi_Integer] :=
  ProjectiveQuadrantSizeOp[m] - Abs[ProjectiveWrapOp[m, phi]];

ProjectiveCharSOp[m_Integer, phi_Integer] :=
  ProjectiveCharCOp[m, phi - ProjectiveQuadrantSizeOp[m]];

(* Branchless Gray-code semidirect butterfly Z_{2^m} =~ Z_4^Gray |>< Z_{2^(m-2)} for any m >= 3 *)
ProjectiveGrayButterflyPairOp[m_Integer, phi_Integer] := Module[
  {nRing = ProjectiveRingSizeOp[m], half = ProjectiveHalfTurnOp[m], qSq = ProjectiveQuadrantSizeOp[m],
   u, q1, q0, rem, comp, cMag, sMag, cNeg, sNeg},
  u    = Mod[phi, nRing];
  q1   = Quotient[u, half];
  q0   = Mod[Quotient[u, qSq], 2];
  rem  = Mod[u, qSq];
  comp = qSq - rem;
  cMag = (1 - q0) * comp + q0 * rem;
  sMag = (1 - q0) * rem  + q0 * comp;
  cNeg = Mod[q1 + q0, 2];
  sNeg = q1;
  {(1 - 2 * cNeg) * cMag, (1 - 2 * sNeg) * sMag}
];

(* Exact discrete Wronskian W_m(phi) === 2^(m-2) for all phi in Z_{2^m} and all m >= 3 *)
ProjectiveWronskianOp[m_Integer, phi_Integer] :=
  ProjectiveCharCOp[m, phi] * ProjectiveCharSOp[m, phi + 1] -
  ProjectiveCharSOp[m, phi] * ProjectiveCharCOp[m, phi + 1];

(* Homogeneous Pythagorean lift [P_x : P_y : P_r] in P^2(Q) on S^1: P_x^2 + P_y^2 === P_r^2 for all m *)
ProjectivePythagoreanTripleOp[m_Integer, phi_Integer] := Module[
  {c = ProjectiveCharCOp[m, phi], s = ProjectiveCharSOp[m, phi],
   q2 = ProjectiveQuadrantSizeOp[m]^2, s2, cNonZero, cPos, sgnC},
  s2 = s^2;
  cNonZero = 1 - Quotient[1, Abs[c] + 1];
  cPos = 1 - Quotient[1, Max[0, c] + 1];
  sgnC = cNonZero * (2 * cPos - 1);
  {sgnC * (q2 - s2), 2 * ProjectiveQuadrantSizeOp[m] * s, q2 + s2}
];

(* Exact rational point (x, y) in S^1(Q) with x^2 + y^2 === 1 and scale-invariance under (m, phi) -> (m+1, 2*phi) *)
ProjectiveUnitCircleRationalPointOp[m_Integer, phi_Integer] := Module[
  {p = ProjectivePythagoreanTripleOp[m, phi]},
  {p[[1]] / p[[3]], p[[2]] / p[[3]]}
];

(* Closed-form N-dimensional quadratic integral over (Z_{2^m})^nDim (dl^nDim) in Z[i] for any m >= 3, nDim >= 1 *)
ProjectiveNDQuadraticExponentOp[m_Integer, nDim_Integer] :=
  (m - 2) + nDim * Quotient[m, 2] + Quotient[nDim, 2] * Mod[m, 2];

ProjectiveNDQuadraticSumClosedOp[m_Integer, nDim_Integer] :=
  ((1 + I) ^ nDim) * (2 ^ ProjectiveNDQuadraticExponentOp[m, nDim]);

Projective1DQuadraticSumDirectOp[m_Integer] := Module[
  {nRing = ProjectiveRingSizeOp[m]},
  Sum[ProjectiveCharCOp[m, x^2] + I * ProjectiveCharSOp[m, x^2], {x, 0, nRing - 1}]
];

(* Exact branch-carry cohomology dphi = 2^m * k_xy + A_xy and contour winding index for any m >= 3 *)
ProjectiveEdgeBranchCarryOp[m_Integer, phiY_Integer, phiX_Integer] := Module[
  {nRing = ProjectiveRingSizeOp[m], half = ProjectiveHalfTurnOp[m], uY, uX, dphi, k, a},
  uY   = Mod[phiY, nRing];
  uX   = Mod[phiX, nRing];
  dphi = uY - uX;
  k    = Floor[(dphi + (half - 1)) / nRing];
  a    = dphi - nRing * k;
  {a, k}
];

ProjectiveContourWindingIndexOp[m_Integer, closedLoopPhases_List] := Module[
  {len = Length[closedLoopPhases]},
  -Sum[
    ProjectiveEdgeBranchCarryOp[m, closedLoopPhases[[Mod[k, len] + 1]], closedLoopPhases[[k]]][[2]],
    {k, 1, len}
  ]
];

(* Bit-exact reversible symplectic Verlet integrator on Z_{2^m} without Heisenberg deadband *)
ProjectiveSymplecticStepForwardOp[m_Integer, phiCurr_Integer, phiPast_Integer, drive_Integer] :=
  {Mod[2 * phiCurr - phiPast + drive, ProjectiveRingSizeOp[m]], Mod[phiCurr, ProjectiveRingSizeOp[m]]};

ProjectiveSymplecticStepReverseOp[m_Integer, phiNext_Integer, phiCurr_Integer, drive_Integer] :=
  {Mod[phiCurr, ProjectiveRingSizeOp[m]], Mod[2 * phiCurr - phiNext + drive, ProjectiveRingSizeOp[m]]};

(* O(d^2) isotropic d-dimensional stencil O_2(d) = {e in {-1,0,1}^d : |e|^2 == 2} for arbitrary d >= 2 *)
ProjectiveIsotropicStencilO2Op[d_Integer] := ProjectiveIsotropicStencilO2Op[d] =
  Select[Tuples[ConstantArray[{-1, 0, 1}, d]], (# . #) == 2 &];

ProjectiveIsotropicSecondMomentOp[d_Integer] := Module[
  {o2 = ProjectiveIsotropicStencilO2Op[d]},
  Table[Sum[o2[[k, i]] * o2[[k, j]], {k, 1, Length[o2]}], {i, 1, d}, {j, 1, d}]
];
