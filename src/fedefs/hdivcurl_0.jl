"""
````
abstract type HDIVCURL0{edim} <: AbstractHdivcurlFiniteElement where {edim<:Int}
````

Hdivcurl-conforming matrix-valued lowest-order element.

allowed ElementGeometries:
- Triangle2D
"""
abstract type HDIVCURL0{edim} <: AbstractHdivcurlFiniteElement where {edim<:Int} end

function Base.show(io::Core.IO, ::Type{<:HDIVCURL0{edim}}) where {edim}
    print(io,"HDIVCURL0{$edim}")
end

get_ncomponents(FEType::Type{<:HDIVCURL0}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}, Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVCURL0}, EG::Type{<:AbstractElementGeometry}) = 1
get_ndofs(::Type{ON_CELLS}, FEType::Type{HDIVCURL0{2}}, EG::Type{<:AbstractElementGeometry}) = num_faces(EG)

get_polynomialorder(::Type{<:HDIVCURL0{2}}, ::Type{<:AbstractElementGeometry1D}) = 0;
get_polynomialorder(::Type{<:HDIVCURL0{2}}, ::Type{<:AbstractElementGeometry2D}) = 1;

get_dofmap_pattern(FEType::Type{<:HDIVCURL0{2}}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "f1"
get_dofmap_pattern(FEType::Type{<:HDIVCURL0{2}}, ::Union{Type{FaceDofs},Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "i1"

isdefined(FEType::Type{<:HDIVCURL0}, ::Type{<:Triangle2D}) = true


# on faces dofs are only continuous tangential-normal fluxes
function get_basis(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}, Type{<:ON_BFACES}, Type{<:ON_FACES}}, ::Type{<:HDIVCURL0}, ::Type{<:AbstractElementGeometry})
    function closure(refbasis, xref)
        refbasis[1,1] = 1
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVCURL0{2}}, ::Type{<:Triangle2D})
    function closure(refbasis, xref)
        refbasis[:,2] = sqrt(2) * [-1, 0, 0, 1]
        refbasis[:,3] = [0.5, 0, 1, -0.5]
        refbasis[:,1] = [0.5, -1, 0, -0.5]
        return nothing
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv,Ti,<:HDIVCURL0,APT}, EG::Type{<:AbstractElementGeometry2D}) where {Tv,Ti,APT}
    xCellFaceSigns = FE.xgrid[CellFaceSigns]
    nfaces = num_faces(EG)
    function closure(coefficients, cell)
        # multiplication with normal vector signs
        for j = 1 : nfaces,  k = 1 : size(coefficients,1)
            coefficients[k,j] = xCellFaceSigns[j,cell];
        end
        return nothing
    end
end   

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv,Ti,<:HDIVCURL0,APT}, EG::Type{<:AbstractElementGeometry3D}) where {Tv,Ti,APT}
    xCellEdgeSigns = FE.xgrid[CellEdgeSigns]
    nedges = num_edges(EG)
    function closure(coefficients, cell)
        # multiplication with normal vector signs
        for j = 1 : nedges,  k = 1 : size(coefficients,1)
            coefficients[k,j] = xCellEdgeSigns[j,cell];
        end
        return nothing
    end
end     
