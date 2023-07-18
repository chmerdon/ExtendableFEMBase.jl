
"""
````
abstract type H1P2H{ncomponents,edim} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int}
````

Continuous piecewise second-order polynomials with hierarchical basis.

allowed ElementGeometries:
- Edge1D (quadratic polynomials)
- Triangle2D (quadratic polynomials)
- Quadrilateral2D (Q2 space)
- Tetrahedron3D (quadratic polynomials)
"""
abstract type H1P2H{ncomponents,edim} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int} end

function Base.show(io::Core.IO, ::Type{<:H1P2H{ncomponents,edim} }) where {ncomponents,edim} 
    print(io,"H1P2H{$ncomponents,$edim}")
end

get_ncomponents(FEType::Type{<:H1P2H}) = FEType.parameters[1]
get_edim(FEType::Type{<:H1P2H}) = FEType.parameters[2]

get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:H1P2H}, EG::Type{<:AbstractElementGeometry0D}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}}, FEType::Type{<:H1P2H}, EG::Type{<:AbstractElementGeometry1D}) = 3*FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1P2H}, EG::Type{<:Union{AbstractElementGeometry1D, Triangle2D, Tetrahedron3D}}) = Int((FEType.parameters[2])*(FEType.parameters[2]+1)/2*FEType.parameters[1])
get_ndofs(::Type{<:ON_CELLS},FEType::Type{<:H1P2H}, EG::Type{<:Union{AbstractElementGeometry1D, Triangle2D, Tetrahedron3D}}) = Int((FEType.parameters[2]+1)*(FEType.parameters[2]+2)/2*FEType.parameters[1])
get_ndofs(::Type{<:ON_CELLS},FEType::Type{<:H1P2H}, EG::Type{<:Quadrilateral2D}) = 8*FEType.parameters[1]

get_polynomialorder(::Type{<:H1P2H}, ::Type{<:Edge1D}) = 2;
get_polynomialorder(::Type{<:H1P2H}, ::Type{<:Triangle2D}) = 2;
get_polynomialorder(::Type{<:H1P2H}, ::Type{<:Quadrilateral2D}) = 3;
get_polynomialorder(::Type{<:H1P2H}, ::Type{<:Tetrahedron3D}) = 2;

get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"
get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "N1F1"
get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry3D}) = "N1E1"
get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Union{Type{FaceDofs},Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry0D}) = "N1"
get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Union{Type{FaceDofs},Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"
get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Union{Type{FaceDofs},Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry2D}) = "N1E1"
get_dofmap_pattern(FEType::Type{<:H1P2H}, ::Type{EdgeDofs}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"

isdefined(FEType::Type{<:H1P2H}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1P2H}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1P2H}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1P2H}, ::Type{<:Tetrahedron3D}) = true

interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1P2H{ncomponents,edim}}, ::Type{Edge1D}) where {ncomponents,edim} = 2

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{AT_NODES}, exact_function!; items = [], bonus_quadorder::Int = 0, time = 0) where {Tv,Ti,FEType <: H1P2H,APT}
    edim = get_edim(FEType)
    nnodes = size(FE.xgrid[Coordinates],2)
    if edim == 1
        nedges = num_sources(FE.xgrid[CellNodes])
    elseif edim == 2
        nedges = num_sources(FE.xgrid[FaceNodes])
    elseif edim == 3
        nedges = num_sources(FE.xgrid[EdgeNodes])
    end

    point_evaluation!(Target, FE, AT_NODES, exact_function!; items = items, component_offset = nnodes + nedges, time = time)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_EDGES}, exact_function!; items = [], bonus_quadorder::Int = 0, time = 0) where {Tv,Ti,FEType <: H1P2H,APT}
    edim = get_edim(FEType)
    if edim == 3
        # delegate edge nodes to node interpolation
        subitems = slice(FE.xgrid[EdgeNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, time = time)

        # perform edge mean interpolation
        ensure_moments!(Target, FE, ON_EDGES, exact_function!; items = items, time = time)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_FACES}, exact_function!; items = [], bonus_quadorder::Int = 0, time = 0) where {Tv,Ti,FEType <: H1P2H,APT}
    edim = get_edim(FEType)
    if edim == 2
        # delegate face nodes to node interpolation
        subitems = slice(FE.xgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, time = time)

        # perform face mean interpolation
        ensure_moments!(Target, FE, ON_FACES, exact_function!; items = items, time = time)
    elseif edim == 3
        # delegate face edges to edge interpolation
        subitems = slice(FE.xgrid[FaceEdges], items)
        interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, time = time)
    elseif edim == 1
        # delegate face nodes to node interpolation
        subitems = slice(FE.xgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, time = time)
    end
end


function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_CELLS}, exact_function!; items = [], bonus_quadorder::Int = 0, time = 0) where {Tv,Ti,FEType <: H1P2H,APT}
    edim = get_edim(FEType)
    if edim == 2
        # delegate cell faces to face interpolation
        subitems = slice(FE.xgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, time = time)
    elseif edim == 3
        # delegate cell edges to edge interpolation
        subitems = slice(FE.xgrid[CellEdges], items)
        interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, time = time)
    elseif edim == 1
        # delegate cell nodes to node interpolation
        subitems = slice(FE.xgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, time = time)

        # preserve cell integral
        ensure_moments!(Target, FE, ON_CELLS, exact_function!; items = items, time = time)
    end
end


function get_basis(::Type{<:AssemblyType},FEType::Type{H1P2H{ncomponents,edim}}, ::Type{<:Vertex0D}) where {ncomponents,edim}
    function closure(refbasis,xref)
        for k = 1 : ncomponents
            refbasis[k,k] = 1
        end
    end
end

function get_basis(::Type{<:AssemblyType},FEType::Type{H1P2H{ncomponents,edim}}, ::Type{<:Edge1D}) where {ncomponents,edim}
    function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1]
        for k = 1 : ncomponents
            refbasis[3*k-2,k] = 1 - xref[1] # node 1
            refbasis[3*k-1,k] = xref[1]     # node 2
            refbasis[3*k,k] = 4*refbasis[end]*xref[1]                       # face 1
        end
    end
end

function get_basis(::Type{<:AssemblyType},FEType::Type{H1P2H{ncomponents,edim}}, ::Type{<:Triangle2D}) where {ncomponents,edim}
    function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1] - xref[2] # store last barycentric coordinate
        for k = 1 : ncomponents
            refbasis[6*k-5,k] = 1 - xref[1] -xref[2]      # node 1
            refbasis[6*k-4,k] = xref[1]                   # node 2
            refbasis[6*k-3,k] = xref[2]                   # node 3
            refbasis[6*k-2,k] = 4*refbasis[end]*xref[1]                     # face 1
            refbasis[6*k-1,k] = 4*xref[1]*xref[2]                           # face 2
            refbasis[6*k,k] = 4*xref[2]*refbasis[end]                       # face 3
        end
    end
end
