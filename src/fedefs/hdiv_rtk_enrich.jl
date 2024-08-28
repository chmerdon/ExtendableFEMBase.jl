"""
````
abstract type HDIVRTkENRICH{k,edim} <: AbstractHdivFiniteElement where {edim<:Int}
````

Internal (normal-zero) Hdiv-conforming vector-valued (ncomponents = edim) Raviart-Thomas space of order k ≥ 1
with the additional orthogonality property that their divergences are L2-orthogonal on P_{k-edim+1}.
Example: HDIVRTkENRICH{1,2} gives the edim interior RT1 bubbles (= normal-trace-free) on a triangle, their divergences
have integral mean zero; HDIVRTkENRICH{2,2} gives three RT2 bubbles on a triangle whose divergences are L2-orthogonal
onto all P1 functions. The maximal order for k is 4 on a Triangle2D (edim = 2) and 3 on Tetrahedron3D (edim = 3).
These spaces have no approximation power on their own, but can be used as enrichment spaces in divergence-free
schemes for incompressible Stokes problems.

allowed ElementGeometries:
- Triangle2D
- Tetrahedron3D
"""
abstract type HDIVRTkENRICH{edim, k, stack} <: AbstractHdivFiniteElement where {k <: Int, edim <: Int, stack <: Bool} end

const _num_RTk_enrich_bubbles = [[2, 3, 4], [3, 9]]
const _num_RTk_enrich_bubbles_stacked = [[2, 5, 9], [3, 9]]
get_ncomponents(::Type{<:HDIVRTkENRICH{edim, order}}) where {order, edim} = edim
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRTkENRICH{edim, order, false}}, EG::Type{<:AbstractElementGeometry}) where {edim, order} = _num_RTk_enrich_bubbles[edim-1][order]
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRTkENRICH{edim, order, true}}, EG::Type{<:AbstractElementGeometry}) where {edim, order} = _num_RTk_enrich_bubbles_stacked[edim-1][order]
get_polynomialorder(::Type{<:HDIVRTkENRICH{edim, order}}, ::Type{<:AbstractElementGeometry}) where {edim, order} = order + 1
get_dofmap_pattern(FEType::Type{<:HDIVRTkENRICH}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "i$(get_ndofs(ON_CELLS, FEType, EG))"

isdefined(FEType::Type{<:HDIVRTkENRICH}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:HDIVRTkENRICH}, ::Type{<:Tetrahedron3D}) = true

## basis on reference triangle
function get_basis(::Type{ON_CELLS}, ::Type{HDIVRTkENRICH{2, order, false}}, EG::Type{<:Triangle2D}) where {order}
	@assert order in 1:3
	#@show chosen_bubbles
	function closure(refbasis, xref)
		# RT1 bubbles
		# ψ_j is RT0 basis function of j-th edge
		# φ_j is nodal basis function
		# node [3,1,2] is opposite to face [1,2,3]

		# 3
		# |\
		# | \ E2
		# |  \
		# |___\
		# 1 E1 2

		# RT1 cell bubbles
		refbasis[1, 1] = xref[2] * xref[1]
		refbasis[1, 2] = xref[2] * (xref[2] - 1)  # = φ_3 ψ_1
		refbasis[2, 1] = xref[1] * (xref[1] - 1)
		refbasis[2, 2] = xref[1] * xref[2]  # = φ_2 ψ_3

		# refbasis[1,1] = 8*(1-xref[1]-xref[2])*xref[1] - 8*xref[1]*(xref[1]-1)
		# refbasis[1,2] = 8*(1-xref[1]-xref[2])*xref[2] - 8*xref[1]*(xref[2])

		#  refbasis[2,1] = 8*xref[1]*(xref[1]-1) - 8*xref[2]*xref[1]
		#  refbasis[2,2] = 8*xref[1]*(xref[2]) - 8*xref[2]*(xref[2]-1)

		# minimal selection with L^2 orthogonality: div(*) ⟂ P_{order-1})
		if order == 1
			# just the RT1 bubbles above
		elseif order == 2
			# special bubbles with zero lowest order divergence moments
			for k ∈ 1:2
				refbasis[3, k] = (5 * (1 - xref[1] - xref[2]) - 2) * (-refbasis[1, k] - refbasis[2, k])
				refbasis[1, k] = (5 * xref[2] - 2) * refbasis[1, k]
				refbasis[2, k] = (5 * xref[1] - 2) * refbasis[2, k]
			end
		elseif order == 3
			for k ∈ 1:2
				refbasis[3, k] = (7 * (1 - xref[1] - xref[2])^2 - 6 * (1 - xref[1] - xref[2]) + 1) * (-refbasis[1, k] - refbasis[2, k]) / 7
				refbasis[4, k] = -2 * xref[1] * xref[2] * refbasis[2, k] + 2 // 45 * (-refbasis[1, k] + 4 * refbasis[2, k])
				refbasis[2, k] = (7 * xref[1]^2 - 6 * xref[1] + 1) * refbasis[2, k] / 7
				refbasis[1, k] = (7 * xref[2]^2 - 6 * xref[2] + 1) * refbasis[1, k] / 7
				refbasis[4, k] += (3 * refbasis[1, k] + 2 * refbasis[2, k] - 3 * refbasis[3, k]) / 70
			end
		end
	end
end


## basis on reference tetrahedron
function get_basis(::Type{ON_CELLS}, ::Type{HDIVRTkENRICH{3, order, false}}, EG::Type{<:Tetrahedron3D}) where {order}
	@assert order in 1:2
	#@show chosen_bubbles
	function closure(refbasis, xref)
		# all RT1 bubbles
		refbasis[1, 1] = 2 * xref[3] * xref[1]
		refbasis[1, 2] = 2 * xref[3] * xref[2]
		refbasis[1, 3] = 2 * xref[3] * (xref[3] - 1)
		refbasis[2, 1] = 2 * xref[2] * xref[1]
		refbasis[2, 2] = 2 * xref[2] * (xref[2] - 1)
		refbasis[2, 3] = 2 * xref[2] * xref[3]
		refbasis[3, 1] = 2 * xref[1] * (xref[1] - 1)
		refbasis[3, 2] = 2 * xref[1] * xref[2]
		refbasis[3, 3] = 2 * xref[1] * xref[3]

		if order == 1
			# nothing to add (but enrichment need additional RT0 handled by seperate FESpace/FEVectorBlock)
		elseif order == 2
			for k ∈ 1:3
				refbasis[4, k] = (6 * (1 - xref[1] - xref[2] - xref[3]) - 1) * refbasis[3, k] + (6 * xref[1] - 1) * (-refbasis[1, k] - refbasis[2, k] - refbasis[3, k]) # (1,2)
				refbasis[5, k] = (6 * (1 - xref[1] - xref[2] - xref[3]) - 1) * refbasis[2, k] + (6 * xref[2] - 1) * (-refbasis[1, k] - refbasis[2, k] - refbasis[3, k]) # (1,3)
				refbasis[6, k] = (6 * (1 - xref[1] - xref[2] - xref[3]) - 1) * refbasis[1, k] + (6 * xref[3] - 1) * (-refbasis[1, k] - refbasis[2, k] - refbasis[3, k]) # (1,4)
				refbasis[7, k] = (6 * xref[1] - 1) * refbasis[2, k] + (6 * xref[2] - 1) * refbasis[3, k] # (2,3)
				refbasis[8, k] = (6 * xref[1] - 1) * refbasis[1, k] + (6 * xref[3] - 1) * refbasis[3, k] # (2,4)
				refbasis[9, k] = (6 * xref[2] - 1) * refbasis[1, k] + (6 * xref[3] - 1) * refbasis[2, k] # (3,4)
			end
		end
	end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRTkENRICH{2, order, true}}, EG::Type{<:Triangle2D}) where {order}
	@assert order in 1:3
	refbasisRT0 = zeros(Real, 3, 2)
	refbasisRT1 = zeros(Real, 2, 2) # only bubbles
	#@show chosen_bubbles
	function closure(refbasis, xref)
		refbasisRT0[1, 1] = xref[1]
		refbasisRT0[1, 2] = xref[2] - 1
		refbasisRT0[2, 1] = xref[1]
		refbasisRT0[2, 2] = xref[2]
		refbasisRT0[3, 1] = xref[1] - 1
		refbasisRT0[3, 2] = xref[2]

		# all RT 1 bubbles
		refbasisRT1[1, :] .= xref[2] * refbasisRT0[1, :]    # = φ_3 ψ_1
		refbasisRT1[2, :] .= xref[1] * refbasisRT0[3, :]    # = φ_2 ψ_3

		# selection
		next = 1
		if order >= 1
			#refbasis[1:total_bubbles,:] .= refbasisRT1[chosen_bubbles,:]
			refbasis[next, :] = refbasisRT1[1, :]
			next += 1
			refbasis[next, :] = refbasisRT1[2, :]
			next += 1
		end
		if order >= 2
			# special bubbles with zero lowest order divergence moments
			refbasis[next, :] .= (5 * (1 - xref[1] - xref[2]) - 2) * (-refbasisRT1[1, :] - refbasisRT1[2, :])
			next += 1
			refbasis[next, :] .= (5 * xref[1] - 2) * refbasisRT1[2, :]
			next += 1
			refbasis[next, :] .= (5 * xref[2] - 2) * refbasisRT1[1, :]
			next += 1
		end
		if order >= 3
			# special bubbles with zero lowest order divergence moments
			refbasis[next, :] .= -(7 * (1 - xref[1] - xref[2])^2 - 6 * (1 - xref[1] - xref[2]) + 1) * (-refbasisRT1[1, :] - refbasisRT1[2, :]) / 7
			next += 1
			refbasis[next, :] .= (7 * xref[1]^2 - 6 * xref[1] + 1) * refbasisRT1[2, :] / 7
			next += 1
			refbasis[next, :] .= (7 * xref[2]^2 - 6 * xref[2] + 1) * refbasisRT1[1, :] / 7
			next += 1
			refbasis[next, :] .= -2 * xref[1] * xref[2] * refbasisRT1[2, :] + 2 // 45 * (-refbasisRT1[1, :] - refbasisRT1[2, :] + 5 * refbasisRT1[2, :]) + (3 * refbasis[3, :] + 2 * refbasis[4, :] - 3 * refbasis[5, :]) / 70
		end
	end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRTkENRICH{3, order, true}}, EG::Type{<:Tetrahedron3D}) where {order}
	@assert order in 1:3
	refbasisRT0 = zeros(Real, 4, 3)
	refbasisRT1 = zeros(Real, 3, 3) # only bubbles
	#@show chosen_bubbles
	function closure(refbasis, xref)
		# RT0
		refbasisRT0[1, 1] = 2 * xref[1]
		refbasisRT0[1, 2] = 2 * xref[2]
		refbasisRT0[1, 3] = 2 * (xref[3] - 1)
		refbasisRT0[2, 1] = 2 * xref[1]
		refbasisRT0[2, 2] = 2 * (xref[2] - 1)
		refbasisRT0[2, 3] = 2 * xref[3]
		refbasisRT0[3, 1] = 2 * xref[1]
		refbasisRT0[3, 2] = 2 * xref[2]
		refbasisRT0[3, 3] = 2 * xref[3]
		refbasisRT0[4, 1] = 2 * (xref[1] - 1)
		refbasisRT0[4, 2] = 2 * xref[2]
		refbasisRT0[4, 3] = 2 * xref[3]

		# all RT1 bubbles
		refbasisRT1[1, :] .= xref[3] * refbasisRT0[1, :] # = φ_4 ψ_1
		refbasisRT1[2, :] .= xref[2] * refbasisRT0[2, :] # = φ_3 ψ_2
		refbasisRT1[3, :] .= xref[1] * refbasisRT0[4, :] # = φ_2 ψ_4

		if order == 1
			refbasis[1:3, :] .= refbasisRT1
		elseif order >= 2
			refbasis[1:3, :] .= refbasisRT1
			refbasis[4, :] .= (6 * (1 - xref[1] - xref[2] - xref[3]) - 1) * refbasisRT1[3, :] + (6 * xref[1] - 1) * (-refbasisRT1[1, :] - refbasisRT1[2, :] - refbasisRT1[3, :]) # (1,2)
			refbasis[5, :] .= (6 * (1 - xref[1] - xref[2] - xref[3]) - 1) * refbasisRT1[2, :] + (6 * xref[2] - 1) * (-refbasisRT1[1, :] - refbasisRT1[2, :] - refbasisRT1[3, :]) # (1,3)
			refbasis[6, :] .= (6 * (1 - xref[1] - xref[2] - xref[3]) - 1) * refbasisRT1[1, :] .+ (6 * xref[3] - 1) * (-refbasisRT1[1, :] - refbasisRT1[2, :] - refbasisRT1[3, :]) # (1,4)
			refbasis[7, :] .= (6 * xref[1] - 1) * refbasisRT1[2, :] + (6 * xref[2] - 1) * refbasisRT1[3, :] # (2,3)
			refbasis[8, :] .= (6 * xref[1] - 1) * refbasisRT1[1, :] + (6 * xref[3] - 1) * refbasisRT1[3, :] # (2,4)
			refbasis[9, :] .= (6 * xref[2] - 1) * refbasisRT1[1, :] + (6 * xref[3] - 1) * refbasisRT1[2, :] # (3,4)
		end
	end
end
