module ExtendableFEMBase

using ExtendableGrids # + some reexports from there
export Edge1D, Triangle2D, Parallelogram2D, Tetrahedron3D, Parallelepiped3D
export ON_CELLS, ON_FACES, ON_IFACES, ON_BFACES, ON_EDGES, ON_BEDGES
using ExtendableSparse
using SparseArrays
using DiffResults
using LinearAlgebra
using ForwardDiff
using DocStringExtensions
using Printf
using UnicodePlots
using Term
using SpecialPolynomials
using Polynomials


include("qpinfos.jl")
export QPInfos

include("quadrature.jl")
export QuadratureRule
export VertexRule
export integrate!, integrate, ref_integrate!

include("functionoperators.jl")
export AbstractFunctionOperator
export StandardFunctionOperator
export Identity, IdentityComponent
export NormalFlux, TangentFlux
export Gradient
export SymmetricGradient, TangentialGradient
export Divergence
export CurlScalar, Curl2D, Curl3D
export Laplacian, Hessian, SymmetricHessian
export Trace, Deviator
export NeededDerivative4Operator, Length4Operator, QuadratureOrderShift4Operator, DefaultName4Operator
export OperatorPair, OperatorTriple


include("finiteelements.jl") # also includes dofmaps.jl and feevaluator*.jl
export DofMap
export CellDofs, FaceDofs, EdgeDofs, BFaceDofs, BEdgeDofs
export CellDofsParent, FaceDofsParent, EdgeDofsParent, BFaceDofsParent, BEdgeDofsParent
export DofMapTypes
export Dofmap4AssemblyType, ItemGeometries4DofMap, EffAT4AssemblyType, ParentDofmap4Dofmap
export AbstractFiniteElement
export FESpace, FESpaces, get_AT, get_FEType
export boundarydofs, ndofs, broken

export AbstractH1FiniteElement
export H1BUBBLE, L2P0, H1P1, H1P2, H1P2B, H1MINI, H1CR, H1P3, H1Pk
export L2P1
export H1Q1, H1Q2

export AbstractH1FiniteElementWithCoefficients
export H1BR, H1P1TEB

export AbstractHdivFiniteElement
export HDIVRT0, HDIVBDM1, HDIVRT1, HDIVBDM2, HDIVRTk
export HDIVRTkENRICH

export AbstractHcurlFiniteElement
export HCURLN0, HCURLN1

export get_AT
export get_polynomialorder, get_ndofs, get_ndofs_all
export get_ncomponents, get_edim
export get_basis, get_coefficients, get_basissubset

export interpolate! # must be defined separately by each FEdefinition
export nodevalues, continuify
export nodevalues!, nodevalues_subset!
export nodevalues_view

export interpolator_matrix

export FEVectorBlock, FEVector
export dot, norm, norms
export FEMatrixBlock, FEMatrix, _addnz
export fill!, addblock!, addblock_matmul!, lrmatmul, mul!, add!, apply_penalties!
export submatrix

export displace_mesh, displace_mesh!


include("plots.jl")
export unicode_gridplot, unicode_scalarplot

include("reconstructionhandlers.jl")
export ReconstructionHandler, get_rcoefficients!

include("feevaluator.jl")
export FEEvaluator, update_basis!, eval_febe!

include("reconstructionoperators.jl")
export Reconstruct

include("accumvector.jl")
export AccumulatingVector


#
# Print default dict for solver parameters into docstrings
#
function _myprint(dict::Dict{Symbol, Tuple{Any, String}})
	lines_out = IOBuffer()
	for (k, v) in dict
		if typeof(v[1]) <: String
			println(lines_out, "  - $(k): $(v[2]). Default: ''$(v[1])''\n")
		else
			println(lines_out, "  - $(k): $(v[2]). Default: $(v[1])\n")
		end
	end
	String(take!(lines_out))
end
#
# Update solver params from dict
#
function _update_params!(parameters, kwargs)
	for (k, v) in kwargs
		parameters[Symbol(k)] = v
	end
	return nothing
end

include("segment_integrator.jl")
export SegmentIntegrator, initialize!, integrate_segment!

include("point_evaluator.jl")
export PointEvaluator, evaluate!, evaluate_bary!, eval_func, eval_func_bary

include("lazy_interpolate.jl")
export lazy_interpolate!

end # module ExtendableFEMBase.
