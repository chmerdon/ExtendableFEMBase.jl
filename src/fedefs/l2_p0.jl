"""
````
abstract type L2P0{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Piecewise constant polynomials on cells.

allowed on every ElementGeometry
"""
abstract type L2P0{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int} end

function Base.show(io::Core.IO, ::Type{<:L2P0{ncomponents}}) where {ncomponents}
    print(io,"L2P0{$ncomponents}")
end

get_ncomponents(FEType::Type{<:L2P0}) = FEType.parameters[1]
get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:L2P0}, EG::Type{<:AbstractElementGeometry}) = FEType.parameters[1]

get_polynomialorder(::Type{<:L2P0}, ::Type{<:AbstractElementGeometry}) = 0

get_dofmap_pattern(FEType::Type{<:L2P0}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "I1"
get_dofmap_pattern(FEType::Type{<:L2P0}, ::Union{Type{FaceDofs},Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "C1"

isdefined(FEType::Type{<:L2P0}, ::Type{<:AbstractElementGeometry}) = true

function ExtendableGrids.interpolate!(Target::AbstractArray{T,1}, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_CELLS}, exact_function!; items = [], time = time) where {T,Tv,Ti,FEType <: L2P0,APT}
    xCellVolumes = FE.xgrid[CellVolumes]
    ncells = num_sources(FE.xgrid[CellNodes])
    if items == []
        items = 1 : ncells
    else
        items = filter(!iszero, items)
    end
    ncomponents = get_ncomponents(FEType)
    integrals4cell = zeros(T,ncomponents,ncells)
    integrate!(integrals4cell, FE.xgrid, ON_CELLS, exact_function!; items = items, time = time)
    for cell in items
        if cell != 0
            for c = 1 : ncomponents
                Target[(cell-1)*ncomponents + c] = integrals4cell[c, cell] / xCellVolumes[cell]
            end
        end
    end    
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_FACES}, exact_function!; items = [], time = 0) where {Tv,Ti,FEType <: L2P0,APT}
    # delegate to node cell interpolation
    subitems = slice(FE.xgrid[FaceCells], items)
    interpolate!(Target, FE, ON_CELLS, exact_function!; items = subitems, time = time)
end

function nodevalues!(Target::AbstractArray{<:Real,2}, Source::AbstractArray{<:Real,1}, FE::FESpace{<:L2P0})
    xCoords = FE.xgrid[Coordinates]
    xCellNodes = FE.xgrid[CellNodes]
    xNodeCells = atranspose(xCellNodes)
    FEType = eltype(FE)
    ncomponents = get_ncomponents(FEType)
    value = 0.0
    nneighbours = 0
    for node = 1 : num_sources(xCoords)
        for c = 1 : ncomponents
            value = 0.0
            nneighbours = num_targets(xNodeCells,node)
            for n = 1 : nneighbours
                value += Source[(xNodeCells[n,node]-1)*ncomponents+c]
            end
            value /= nneighbours
            Target[c,node] = value
        end
    end    
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{L2P0{ncomponents}}, ::Type{<:AbstractElementGeometry}) where {ncomponents}
    function closure(refbasis, xref)
        for k = 1 : ncomponents
            refbasis[k,k] = 1.0
        end
    end
end