(* ::Package:: *)
(* Root entry point: loads pure M-layer Core CA (Core/*.wl) and T<->M Lab layers above it (Lab/*.wl) *)

mtCarrierBaseDir = DirectoryName[$InputFileName];
mtCoreDir        = FileNameJoin[{mtCarrierBaseDir, "Core"}];
mtLabDir         = FileNameJoin[{mtCarrierBaseDir, "Lab"}];

(* Layer 1: Pure M-Layer Cellular Automaton Core (0 domain/SI/biology/sequence knowledge) *)
Get[FileNameJoin[{mtCoreDir, "Spec.wl"}]];
Get[FileNameJoin[{mtCoreDir, "FaceAlgebra.wl"}]];
Get[FileNameJoin[{mtCoreDir, "HolonomyMaxwell.wl"}]];
Get[FileNameJoin[{mtCoreDir, "CarrierAutomaton3D.wl"}]];

(* Layer 2: T <-> M Lab & Domain Adapters above Core *)
Get[FileNameJoin[{mtLabDir, "MToT.wl"}]];
Get[FileNameJoin[{mtLabDir, "ContinuumPhysics.wl"}]];
Get[FileNameJoin[{mtLabDir, "BioCosmology.wl"}]];
Get[FileNameJoin[{mtLabDir, "NuclearShells.wl"}]];
Get[FileNameJoin[{mtLabDir, "ParticleMaterials.wl"}]];
Get[FileNameJoin[{mtLabDir, "ReversibleCalculus.wl"}]];
Get[FileNameJoin[{mtLabDir, "HardwarePhotovoltaics.wl"}]];
Get[FileNameJoin[{mtLabDir, "VerifySuite.wl"}]];

dSolutions = SolveUniqueCatalanDimensionOp[8];
dSolved    = First[dSolutions];
