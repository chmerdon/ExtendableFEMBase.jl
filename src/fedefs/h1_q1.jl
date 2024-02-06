
"""
````
abstract type Q1P1{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Continuous piecewise first-order polynomials on simplices and quads, can be used for mixed geometries.

allowed ElementGeometries:
- Edge1D (P1 space)
- Triangle2D (P1 space)
- Quadrilateral2D (Q1 space)
- Tetrahedron3D (P1 space)
- Hexahedron3D (Q1 space)
"""
abstract type H1Q1{ncomponents} <: AbstractH1FiniteElement where {ncomponents <: Int} end


function Base.show(io::Core.IO, ::Type{<:H1Q1{ncomponents}}) where {ncomponents}
	print(io, "H1Q1{$ncomponents}")
end

get_ncomponents(FEType::Type{<:H1Q1}) = FEType.parameters[1] # is this okay?
get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:H1Q1}, EG::Type{<:AbstractElementGeometry}) = num_nodes(EG) * FEType.parameters[1]

get_polynomialorder(::Type{<:H1Q1}, ::Type{<:Edge1D}) = 1;
get_polynomialorder(::Type{<:H1Q1}, ::Type{<:Triangle2D}) = 1;
get_polynomialorder(::Type{<:H1Q1}, ::Type{<:Tetrahedron3D}) = 1;
get_polynomialorder(::Type{<:H1Q1}, ::Type{<:Quadrilateral2D}) = 2;
get_polynomialorder(::Type{<:H1Q1}, ::Type{<:Hexahedron3D}) = 3;

get_dofmap_pattern(FEType::Type{<:H1Q1}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "N1"
get_dofmap_pattern(FEType::Type{<:H1Q1}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "N1"

isdefined(FEType::Type{<:H1Q1}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1Q1}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1Q1}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1Q1}, ::Type{<:Tetrahedron3D}) = true
isdefined(FEType::Type{<:H1Q1}, ::Type{<:Hexahedron3D}) = true

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q1, APT}
	nnodes = size(FE.dofgrid[Coordinates], 2)
	point_evaluation!(Target, FE, AT_NODES, exact_function; items = items, component_offset = nnodes, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q1, APT}
	# delegate edge nodes to node interpolation
	subitems = slice(FE.dofgrid[EdgeNodes], items)
	interpolate!(Target, FE, AT_NODES, exact_function; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q1, APT}
	# delegate face nodes to node interpolation
	subitems = slice(FE.dofgrid[FaceNodes], items)
	interpolate!(Target, FE, AT_NODES, exact_function; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q1, APT}
	if FE.broken == true
		# broken interpolation
		point_evaluation_broken!(Target, FE, ON_CELLS, exact_function; items = items, kwargs...)
	else
		# delegate cell nodes to node interpolation
		subitems = slice(FE.dofgrid[CellNodes], items)
		interpolate!(Target, FE, AT_NODES, exact_function; items = subitems, kwargs...)
	end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q1{ncomponents}}, ET::Type{<:Union{Vertex0D, AbstractElementGeometry1D, Triangle2D, Tetrahedron3D}}) where {ncomponents}
	edim::Int = dim_element(ET)
	function closure(refbasis, xref)
		for k ∈ 1:ncomponents
			refbasis[(edim+1)*k-edim, k] = 1
			for j ∈ 1:edim
				refbasis[(edim+1)*k-edim, k] -= xref[j]
				refbasis[(edim+1)*k-edim+j, k] = xref[j]
			end
		end
	end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q1{ncomponents}}, ::Type{<:Quadrilateral2D}) where {ncomponents}
	function closure(refbasis, xref)
		refbasis[1, 1] = 1 - xref[1]
		refbasis[2, 1] = 1 - xref[2]

		refbasis[3, 1] = xref[1] * xref[2]
		refbasis[4, 1] = xref[2] * refbasis[1, 1]
		refbasis[1, 1] = refbasis[1, 1] * refbasis[2, 1]
		refbasis[2, 1] = xref[1] * refbasis[2, 1]

		for k ∈ 2:ncomponents, j ∈ 1:4
			refbasis[4*(k-1)+j, k] = refbasis[j, 1]
		end
	end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q1{ncomponents}}, ::Type{<:Hexahedron3D}) where {ncomponents}
	function closure(refbasis, xref)
		refbasis[1, 1] = 1 - xref[1]
		refbasis[2, 1] = 1 - xref[2]
		refbasis[3, 1] = 1 - xref[3]
		refbasis[4, 1] = xref[2] * refbasis[1, 1] * refbasis[3, 1]
		refbasis[5, 1] = xref[3] * refbasis[1, 1] * refbasis[2, 1]
		refbasis[7, 1] = xref[1] * xref[2] * xref[3]
		refbasis[6, 1] = xref[1] * refbasis[2, 1] * xref[3]
		refbasis[8, 1] = refbasis[1, 1] * xref[2] * xref[3]
		refbasis[1, 1] = refbasis[1, 1] * refbasis[2, 1] * refbasis[3, 1]
		refbasis[2, 1] = xref[1] * refbasis[2, 1] * refbasis[3, 1]
		refbasis[3, 1] = xref[1] * xref[2] * refbasis[3, 1]
		for k ∈ 2:ncomponents, j ∈ 1:8
			refbasis[8*(k-1)+j, k] = refbasis[j, 1]
		end
	end
end
