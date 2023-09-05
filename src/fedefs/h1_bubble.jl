"""
````
abstract type H1BUBBLE{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Piecewise bubbles (=zero at boundary)

allowed element geometries:
- Edge1D (one quadratic bubble)
- Triangle2D (one cubic bubble)
- Quadrilateral2D (one quartic bubble)
- Tetrahedron3D (one cubic bubble)
"""
abstract type H1BUBBLE{ncomponents} <: AbstractH1FiniteElement where {ncomponents <: Int} end

function Base.show(io::Core.IO, ::Type{<:H1BUBBLE{ncomponents}}) where {ncomponents}
	print(io, "H1BUBBLE{$ncomponents}")
end

get_ncomponents(FEType::Type{<:H1BUBBLE}) = FEType.parameters[1]
get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:H1BUBBLE}, EG::Type{<:AbstractElementGeometry}) = FEType.parameters[1]

get_polynomialorder(::Type{<:H1BUBBLE}, ::Type{<:AbstractElementGeometry1D}) = 2;
get_polynomialorder(::Type{<:H1BUBBLE}, ::Type{<:Triangle2D}) = 3;
get_polynomialorder(::Type{<:H1BUBBLE}, ::Type{<:Quadrilateral2D}) = 4;
get_polynomialorder(::Type{<:H1BUBBLE}, ::Type{<:Tetrahedron3D}) = 4;

get_dofmap_pattern(FEType::Type{<:H1BUBBLE}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "I1"

isdefined(FEType::Type{<:H1BUBBLE}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1BUBBLE}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1BUBBLE}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1BUBBLE}, ::Type{<:Tetrahedron3D}) = true

interior_dofs_offset(::Type{ON_CELLS}, ::Type{H1BUBBLE{ncomponents}}, EG::Type{<:AbstractElementGeometry}) where {ncomponents} = 0

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: H1BUBBLE, APT}
	xCellVolumes = FE.xgrid[CellVolumes]
	ncells = num_sources(FE.xgrid[CellNodes])
	if items == []
		items = 1:ncells
	else
		items = filter(!iszero, items)
	end
	ncomponents = get_ncomponents(FEType)
	integrals4cell = zeros(T, ncomponents, ncells)
	integrate!(integrals4cell, FE.xgrid, ON_CELLS, exact_function!; items = items, kwargs...)
	for cell in items
		if cell != 0
			for c ∈ 1:ncomponents
				Target[(cell-1)*ncomponents+c] = integrals4cell[c, cell] / xCellVolumes[cell]
			end
		end
	end
end

function get_basis(::Type{<:AssemblyType}, ::Type{H1BUBBLE{ncomponents}}, ::Type{<:AbstractElementGeometry1D}) where {ncomponents}
	function closure(refbasis, xref)
		for k ∈ 1:ncomponents
			refbasis[k, k] = 6 * xref[1] * (1 - xref[1])
		end
	end
end

function get_basis(::Type{<:AssemblyType}, ::Type{H1BUBBLE{ncomponents}}, ::Type{<:Triangle2D}) where {ncomponents}
	function closure(refbasis, xref)
		for k ∈ 1:ncomponents
			refbasis[k, k] = 60 * (1 - xref[1] - xref[2]) * xref[1] * xref[2]
		end
	end
end

function get_basis(::Type{<:AssemblyType}, ::Type{H1BUBBLE{ncomponents}}, ::Type{<:Quadrilateral2D}) where {ncomponents}
	function closure(refbasis, xref)
		for k ∈ 1:ncomponents
			refbasis[k, k] = 36 * (1 - xref[1]) * (1 - xref[2]) * xref[1] * xref[2]
		end
	end
end

function get_basis(::Type{<:AssemblyType}, ::Type{H1BUBBLE{ncomponents}}, ::Type{<:Tetrahedron3D}) where {ncomponents}
	function closure(refbasis, xref)
		for k ∈ 1:ncomponents
			refbasis[k, k] = 840 * (1 - xref[1] - xref[2] - xref[3]) * xref[1] * xref[2] * xref[3]
		end
	end
end
