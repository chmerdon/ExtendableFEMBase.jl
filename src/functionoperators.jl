"""
$(TYPEDEF)

root type for FunctionOperators.
"""
abstract type AbstractFunctionOperator end # to dispatch which evaluator of the FE_basis_caller is used

abstract type DiscontinuousFunctionOperator <: AbstractFunctionOperator end
abstract type StandardFunctionOperator <: AbstractFunctionOperator end

abstract type Jump{O} <: DiscontinuousFunctionOperator where {O<:StandardFunctionOperator} end # calculates the jump between both sides of the face
abstract type Average{O} <: DiscontinuousFunctionOperator where {O<:StandardFunctionOperator} end # calculates the average between both sides of the face
abstract type Left{O} <: DiscontinuousFunctionOperator where {O<:StandardFunctionOperator} end # calculates the value on left side of the face
abstract type Right{O} <: DiscontinuousFunctionOperator where {O<:StandardFunctionOperator} end # calculates the value on right side of the face


StandardFunctionOperator(::Type{Jump{O}}) where {O} = O
StandardFunctionOperator(::Type{Average{O}}) where {O} = O
StandardFunctionOperator(::Type{Left{O}}) where {O} = O
StandardFunctionOperator(::Type{Right{O}}) where {O} = O
coeffs(::Type{<:Jump}) = [1,-1]
coeffs(::Type{<:Average}) = [0.5,0.5]
coeffs(::Type{<:Left}) = [1,0]
coeffs(::Type{<:Right}) = [0,1]
is_discontinuous(::Type{<:StandardFunctionOperator}) = false
is_discontinuous(::Type{<:DiscontinuousFunctionOperator}) = true
"""
$(TYPEDEF)

identity operator: evaluates finite element function.
"""
abstract type Identity <: StandardFunctionOperator end # 1*v_h
"""
$(TYPEDEF)

identity operator: evaluates only the c-th component of the finite element function.
"""
abstract type IdentityComponent{c} <: StandardFunctionOperator where {c<:Int} end # 1*v_h(c)
"""
$(TYPEDEF)

identity jump operator: evaluates face jumps of finite element function
"""
#abstract type IdentityDisc{DT<:DiscontinuityTreatment} end # [[v_h]]
"""
$(TYPEDEF)

reconstruction identity operator: evaluates a reconstructed version of the finite element function.

FEreconst specifies the reconstruction space and reconstruction algorithm if it is defined for the finite element that it is applied to.
"""
abstract type AbstractFiniteElement end
abstract type ReconstructionIdentity{FEreconst<:AbstractFiniteElement} <: StandardFunctionOperator end # 1*R(v_h)
#abstract type ReconstructionIdentityDisc{FEreconst<:AbstractFiniteElement, DT<:DiscontinuityTreatment} <: ReconstructionIdentity{FEreconst} end # 1*R(v_h)
"""
$(TYPEDEF)

evaluates the normal-flux of the finite element function.

only available on FACES/BFACES and currently only for H1 and Hdiv elements
"""
abstract type NormalFlux <: StandardFunctionOperator end # v_h * n_F # only for Hdiv/H1 on Faces/BFaceFaces
#abstract type NormalFluxDisc{DT<:DiscontinuityTreatment} <: NormalFlux end # v_h * n_F # only for Hdiv/H1 on Faces/BFaceFaces
"""
$(TYPEDEF)

reconstruction normal flux: evaluates the normal flux of a reconstructed version of the finite element function.

FEreconst specifies the reconstruction space and reconstruction algorithm if it is defined for the finite element that it is applied to.
"""
abstract type ReconstructionNormalFlux{FEreconst<:AbstractFiniteElement} <: StandardFunctionOperator end # R(v_h) * n_F

"""
$(TYPEDEF)

evaluates the tangent-flux of the finite element function.
"""
abstract type TangentFlux <: StandardFunctionOperator end # v_h * t_F # only for HCurlScalar on Edges
"""
$(TYPEDEF)

evaluates the gradient of the finite element function.
"""
abstract type Gradient <: StandardFunctionOperator end # 1*v_h
"""
$(TYPEDEF)

reconstruction gradient operator: evaluates the gradient of a reconstructed version of the finite element function.

FEreconst specifies the reconstruction space and reconstruction algorithm if it is defined for the finite element that it is applied to.
"""
abstract type ReconstructionGradient{FEreconst<:AbstractFiniteElement} <: StandardFunctionOperator end # 1*R(v_h)
"""
$(TYPEDEF)

evaluates the symmetric part of the gradient of the finite element function and returns its Voigt compression
(off-diagonals on position [j,k] and [k,j] are added together and weighted with offdiagval).
"""
abstract type SymmetricGradient{offdiagval} <: StandardFunctionOperator where {offdiagval} end # sym(D_geom(v_h))
"""
$(TYPEDEF)

evaluates the gradient of the tangential part of some vector-valued finite element function.
"""
abstract type TangentialGradient <: StandardFunctionOperator end # D_geom(v_h x n) = only gradient of tangential part of vector-valued function
"""
$(TYPEDEF)

evaluates the Laplacian of some (possibly vector-valued) finite element function.
"""
abstract type Laplacian <: StandardFunctionOperator end # L_geom(v_h)
"""
$(TYPEDEF)

evaluates the full Hessian of some (possibly vector-valued) finite element function.
"""
abstract type Hessian <: StandardFunctionOperator end # D^2(v_h)
"""
$(TYPEDEF)

evaluates only the symmetric part of the Hessian of some (possibly vector-valued) finite element function.
A concatenation of Voigt compressed Hessians for each component of the finite element functions is returned.
The weighting of the mixed derivatives can be steered with offdiagval (e.g. √2 or 1 depending on the use case).
"""
abstract type SymmetricHessian{offdiagval} <: StandardFunctionOperator where {offdiagval} end # sym(D^2(v_h))
"""
$(TYPEDEF)

evaluates the curl of some scalar function in 2D, i.e. the rotated gradient.
"""
abstract type CurlScalar <: StandardFunctionOperator end # only 2D: CurlScalar(v_h) = D(v_h)^\perp

"""
$(TYPEDEF)

evaluates the curl of some two-dimensional vector field, i.e. Curl2D((u1,u2)) = du2/dx1 - du1/dx2
"""
abstract type Curl2D <: StandardFunctionOperator end

"""
$(TYPEDEF)

evaluates the curl of some three-dimensional vector field, i.e. Curl3D(u) = \nabla \times u
"""
abstract type Curl3D <: StandardFunctionOperator end # only 3D: Curl3D(v_h) = D \times v_h
"""
$(TYPEDEF)

evaluates the divergence of the finite element function.
"""
abstract type Divergence <: StandardFunctionOperator end # div(v_h)
"""
$(TYPEDEF)

evaluates the divergence of the reconstructed finite element function.

FEreconst specifies the reconstruction space and reconstruction algorithm if it is defined for the finite element that it is applied to.
"""
abstract type ReconstructionDivergence{FEreconst<:AbstractFiniteElement} <: StandardFunctionOperator end # 1*R(v_h)

abstract type Trace <: StandardFunctionOperator end # tr(v_h)
abstract type Deviator <: StandardFunctionOperator end # dev(v_h)


NeededDerivative4Operator(::Type{<:Identity}) = 0
NeededDerivative4Operator(::Type{<:ReconstructionIdentity}) = 0
NeededDerivative4Operator(::Type{<:IdentityComponent}) = 0
NeededDerivative4Operator(::Type{<:ReconstructionNormalFlux}) = 0
NeededDerivative4Operator(::Type{<:NormalFlux}) = 0
NeededDerivative4Operator(::Type{<:TangentFlux}) = 0
NeededDerivative4Operator(::Type{<:Gradient}) = 1
NeededDerivative4Operator(::Type{<:ReconstructionGradient}) = 1
NeededDerivative4Operator(::Type{<:SymmetricGradient}) = 1
NeededDerivative4Operator(::Type{TangentialGradient}) = 1
NeededDerivative4Operator(::Type{<:Laplacian}) = 2
NeededDerivative4Operator(::Type{<:Hessian}) = 2
NeededDerivative4Operator(::Type{<:SymmetricHessian}) = 2
NeededDerivative4Operator(::Type{CurlScalar}) = 1
NeededDerivative4Operator(::Type{Curl2D}) = 1
NeededDerivative4Operator(::Type{Curl3D}) = 1
NeededDerivative4Operator(::Type{<:Divergence}) = 1
NeededDerivative4Operator(::Type{<:ReconstructionDivergence}) = 1
NeededDerivative4Operator(::Type{Trace}) = 0
NeededDerivative4Operator(::Type{Deviator}) = 0

DefaultName4Operator(::Type{<:StandardFunctionOperator}) = "??"
DefaultName4Operator(::Type{Identity}) = "id"
DefaultName4Operator(::Type{<:ReconstructionIdentity}) = "R"
DefaultName4Operator(IC::Type{<:IdentityComponent{T}}) where T = "id_$(T)"
DefaultName4Operator(::Type{NormalFlux}) = "NormalFlux"
DefaultName4Operator(::Type{<:ReconstructionNormalFlux}) = "R NormalFlux"
DefaultName4Operator(::Type{TangentFlux}) = "TangentialFlux"
DefaultName4Operator(::Type{<:Gradient}) = "∇"
DefaultName4Operator(::Type{<:SymmetricGradient}) = "ϵ"
DefaultName4Operator(::Type{TangentialGradient}) = "TangentialGradient"
DefaultName4Operator(::Type{<:Laplacian}) = "Δ"
DefaultName4Operator(::Type{<:Hessian}) = "H"
DefaultName4Operator(::Type{<:SymmetricHessian}) = "symH"
DefaultName4Operator(::Type{CurlScalar}) = "curl"
DefaultName4Operator(::Type{Curl2D}) = "Curl"
DefaultName4Operator(::Type{Curl3D}) = "∇×"
DefaultName4Operator(::Type{Divergence}) = "div"
DefaultName4Operator(::Type{<:ReconstructionDivergence}) = "div R"
DefaultName4Operator(::Type{<:ReconstructionGradient}) = "∇ R"
DefaultName4Operator(::Type{Trace}) = "tr"
DefaultName4Operator(::Type{Deviator}) = "dev"

function Base.show(io::Core.IO, FO::Type{<:AbstractFunctionOperator})
    print(io,"$(DefaultName4Operator(FO))")
end

# length for operator result
Length4Operator(::Type{<:Identity}, xdim::Int, ncomponents::Int) = ncomponents
Length4Operator(::Type{<:ReconstructionIdentity}, xdim::Int, ncomponents::Int) = ncomponents
Length4Operator(::Type{<:IdentityComponent}, xdim::Int, ncomponents::Int) = 1
Length4Operator(::Type{<:NormalFlux}, xdim::Int, ncomponents::Int) = 1
Length4Operator(::Type{<:ReconstructionNormalFlux}, xdim::Int, ncomponents::Int) = 1
Length4Operator(::Type{<:TangentFlux}, xdim::Int, ncomponents::Int) = 1
Length4Operator(::Type{<:Divergence}, xdim::Int, ncomponents::Int) = Int(ceil(ncomponents/xdim))
Length4Operator(::Type{<:ReconstructionDivergence}, xdim::Int, ncomponents::Int) = Int(ceil(ncomponents/xdim))
Length4Operator(::Type{Trace}, xdim::Int, ncomponents::Int) = ceil(sqrt(ncomponents))
Length4Operator(::Type{CurlScalar}, xdim::Int, ncomponents::Int) = ((xdim == 2) ? xdim*ncomponents : Int(ceil(xdim*(ncomponents/xdim))))
Length4Operator(::Type{Curl2D}, xdim::Int, ncomponents::Int) = 1
Length4Operator(::Type{Curl3D}, xdim::Int, ncomponents::Int) = 3
Length4Operator(::Type{<:Gradient}, xdim::Int, ncomponents::Int) = xdim*ncomponents
Length4Operator(::Type{<:ReconstructionGradient}, xdim::Int, ncomponents::Int) = xdim*ncomponents
Length4Operator(::Type{TangentialGradient}, xdim::Int, ncomponents::Int) = 1
Length4Operator(::Type{<:SymmetricGradient}, xdim::Int, ncomponents::Int) = ((xdim == 2) ? 3 : 6)*Int(ceil(ncomponents/xdim))
Length4Operator(::Type{<:Hessian}, xdim::Int, ncomponents::Int) = xdim*xdim*ncomponents
Length4Operator(::Type{<:SymmetricHessian}, xdim::Int, ncomponents::Int) = ((xdim == 2) ? 3 : 6)*ncomponents
Length4Operator(::Type{<:Laplacian}, xdim::Int, ncomponents::Int) = ncomponents


NeededDerivative4Operator(::Type{Jump{O}}) where {O} = NeededDerivative4Operator(O)
Length4Operator(::Type{Jump{O}}, xdim, nc) where {O} = Length4Operator(O, xdim, nc)
DefaultName4Operator(::Type{Jump{O}}) where {O} = "[[" * DefaultName4Operator(O) * "]]"
NeededDerivative4Operator(::Type{Average{O}}) where {O} = NeededDerivative4Operator(O)
Length4Operator(::Type{Average{O}}, xdim, nc) where {O} = Length4Operator(O, xdim, nc)
DefaultName4Operator(::Type{Average{O}}) where {O} = "{{" * DefaultName4Operator(O) * "}}"
NeededDerivative4Operator(::Type{Left{O}}) where {O} = NeededDerivative4Operator(O)
Length4Operator(::Type{Left{O}}, xdim, nc) where {O} = Length4Operator(O, xdim, nc)
DefaultName4Operator(::Type{Left{O}}) where {O} = DefaultName4Operator(O) * "|_Left"
NeededDerivative4Operator(::Type{Right{O}}) where {O} = NeededDerivative4Operator(O)
Length4Operator(::Type{Right{O}}, xdim, nc) where {O} = Length4Operator(O, xdim, nc)
DefaultName4Operator(::Type{Right{O}}) where {O} = DefaultName4Operator(O) * "|_Right"

QuadratureOrderShift4Operator(operator) = -NeededDerivative4Operator(operator)


# some infrastructure to allow operator pairs (and possibly triplets etc. by recursive application)

"""
````
abstract type OperatorPair{<:StandardFunctionOperator,<:StandardFunctionOperator} <: StandardFunctionOperator
````

allows to evaluate two operators in place of one, e.g. OperatorPair{Identity,Gradient}.
"""
abstract type OperatorPair{OP1 <:StandardFunctionOperator,OP2 <: StandardFunctionOperator} <: StandardFunctionOperator end
Length4Operator(OP::Type{<:OperatorPair}, xdim::Int, ncomponents::Int) = Length4Operator(OP.parameters[1], xdim, ncomponents) + Length4Operator(OP.parameters[2], xdim, ncomponents)
QuadratureOrderShift4Operator(OP::Type{<:OperatorPair}) = max(QuadratureOrderShift4Operator(OP.parameters[1]),QuadratureOrderShift4Operator(OP.parameters[2]))
DefaultName4Operator(OP::Type{<:OperatorPair}) = "(" * DefaultName4Operator(OP.parameters[1]) * "," *  DefaultName4Operator(OP.parameters[2]) * ")"


"""
````
abstract type OperatorTriple{<:StandardFunctionOperator,<:StandardFunctionOperator} <: StandardFunctionOperator
````

allows to evaluate three operators in place of one, e.g. OperatorTriple{Identity,Gradient,Hessian}.
"""
abstract type OperatorTriple{OP1 <:StandardFunctionOperator,OP2 <: StandardFunctionOperator,OP3 <: StandardFunctionOperator} <: StandardFunctionOperator end
Length4Operator(OP::Type{<:OperatorTriple}, xdim::Int, ncomponents::Int) = Length4Operator(OP.parameters[1], xdim, ncomponents) + Length4Operator(OP.parameters[2], xdim, ncomponents) +  + Length4Operator(OP.parameters[3], xdim, ncomponents)
QuadratureOrderShift4Operator(OP::Type{<:OperatorTriple}) = max(max(QuadratureOrderShift4Operator(OP.parameters[1]),QuadratureOrderShift4Operator(OP.parameters[2])),QuadratureOrderShift4Operator(OP.parameters[3]))
DefaultName4Operator(OP::Type{<:OperatorTriple}) = "(" * DefaultName4Operator(OP.parameters[1]) * "," *  DefaultName4Operator(OP.parameters[2]) * "," *  DefaultName4Operator(OP.parameters[3]) * ")"