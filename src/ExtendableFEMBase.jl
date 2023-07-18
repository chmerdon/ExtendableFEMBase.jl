module ExtendableFEMBase

using ExtendableGrids # + some exports from there
using ExtendableSparse
using SuiteSparse
using SparseArrays
export Edge1D, Triangle2D, Parallelogram2D, Tetrahedron3D, Parallelepiped3D
using DiffResults
using LinearAlgebra
using ForwardDiff
using DocStringExtensions
using Printf
using UnicodePlots


include("qpinfos.jl")
export QPInfos

include("quadrature.jl")
export QuadratureRule
export VertexRule
export integrate!, integrate, ref_integrate!

include("functionoperators.jl")
export AbstractFunctionOperator
export DiscontinuousFunctionOperator, StandardFunctionOperator
export Jump, Average, Left, Right, is_discontinuous
export Identity, IdentityComponent, IdentityDisc
export ReconstructionIdentity, ReconstructionIdentityDisc
export ReconstructionGradient, ReconstructionGradientDisc
export ReconstructionDivergence
export ReconstructionNormalFlux
export NormalFlux, NormalFluxDisc, TangentFlux, TangentFluxDisc
export Gradient, GradientDisc
export SymmetricGradient, TangentialGradient
export Divergence, ReconstructionDivergence
export CurlScalar, Curl2D, Curl3D
export Laplacian, Hessian, SymmetricHessian
export Trace, Deviator
export NeededDerivatives4Operator, QuadratureOrderShift4Operator
export Dofmap4AssemblyType, DofitemAT4Operator
export DefaultDirichletBoundaryOperator4FE
export OperatorPair, OperatorTriple
export Δ, ∇, H


include("finiteelements.jl")
export DofMap, CellDofs, FaceDofs, EdgeDofs, BFaceDofs, BEdgeDofs
export DofMapTypes
export AbstractFiniteElement
export FESpace, FESpaces, get_AT, get_FEType
export get_periodic_coupling_info

export AbstractH1FiniteElement
export H1BUBBLE, L2P0, H1P1, H1P2, H1P2B, H1MINI, H1CR, H1P3, H1Pk
export L2P1
export H1Q1, H1Q2

export AbstractH1FiniteElementWithCoefficients
export H1BR, H1P1TEB

export AbstractHdivFiniteElement
export HDIVRT0, HDIVBDM1, HDIVRT1, HDIVRT1INT, HDIVBDM2
export HDIVRTkENRICH

export AbstractHcurlFiniteElement
export HCURLN0

export get_assemblytype
export get_polynomialorder, get_ndofs, get_ndofs_all
export get_ncomponents, get_edim
export get_basis, get_coefficients, get_basissubset
export reconstruct!

export interpolate! # must be defined separately by each FEdefinition
export nodevalues, continuify
export nodevalues!, nodevalues_subset!
export nodevalues_view

export interpolator_matrix

export FEVectorBlock, FEVector
export dot, norm, norms
export FEMatrixBlock, FEMatrix, _addnz
export fill!, addblock!, addblock_matmul!, lrmatmul, mul!, add!, apply_penalties!
export show_entries

export displace_mesh, displace_mesh!

include("reconstructions.jl")
export ReconstructionHandler, get_rcoefficients!

include("feevaluator.jl")
export FEEvaluator, update_basis!, eval_febe!

include("accumvector.jl")
export AccumulatingVector


end # module ExtendableFEMBase.
