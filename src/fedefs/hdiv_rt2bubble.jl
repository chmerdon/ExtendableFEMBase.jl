"""
````
abstract type HDIVRT2BUBBLE{edim,all} <: AbstractHdivFiniteElement where {edim<:Int, all <: Bool}
````

Just internal (normal-zero) Hdiv-conforming vector-valued (ncomponents = edim) Raviart-Thomas space of order 2.
If all == false, only the higher order bubbles are added (in 2D, 4 instead of 6)

allowed ElementGeometries:
- Triangle2D
"""
abstract type HDIVRT2BUBBLE{edim,all} <: AbstractHdivFiniteElement where {edim<:Int,all<:Bool} end

function Base.show(io::Core.IO, ::Type{<:HDIVRT2BUBBLE{edim,all}}) where {edim,all}
    print(io,"HDIVRT2BUBBLE{$edim,$all}")
end

get_ncomponents(FEType::Type{<:HDIVRT2BUBBLE}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRT2BUBBLE}, EG::Type{<:AbstractElementGeometry1D}) = 0
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRT2BUBBLE{2,all}}, EG::Type{<:Triangle2D}) where {all} = all ? 6 : 3

get_polynomialorder(::Type{<:HDIVRT2BUBBLE{2}}, ::Type{<:AbstractElementGeometry1D}) = 2;
get_polynomialorder(::Type{<:HDIVRT2BUBBLE{2}}, ::Type{<:AbstractElementGeometry2D}) = 3;

get_dofmap_pattern(::Type{<:HDIVRT2BUBBLE{2,all}}, ::Type{CellDofs}, EG::Type{<:Triangle2D}) where {all} = all ? "i6" : "i3"

isdefined(::Type{<:HDIVRT2BUBBLE}, ::Type{<:Triangle2D}) = true

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT2BUBBLE{2,all}}, ::Type{<:Triangle2D}) where {all}
    function closure(refbasis,xref)
        refbasis[1,1] = 12*xref[2] * xref[1];        refbasis[1,2] = 12*xref[2] * (xref[2]-1)
        refbasis[2,1] = 12*xref[1] * (xref[1]-1);    refbasis[2,2] = 12*xref[1] * xref[2]
        refbasis[3,1] = 12*xref[2] * xref[1];        refbasis[3,2] = 12*xref[2] * (xref[2]-1)

        refbasis[1,:] .*= xref[1]
        refbasis[2,:] .*= xref[1]
        refbasis[3,:] .*= xref[2]

        if all
            refbasis[4,1] = 12*xref[1] * (xref[1]-1);    refbasis[4,2] = 12*xref[1] * xref[2]
            refbasis[4,:] .*= xref[2]
            refbasis[5,1] = 12*xref[2] * xref[1];        refbasis[5,2] = 12*xref[2] * (xref[2]-1)
            refbasis[6,1] = 12*xref[1] * (xref[1]-1);    refbasis[6,2] = 12*xref[1] * xref[2]
            refbasis[5,:] .*= (1- xref[1] -xref[2])
            refbasis[6,:] .*= (1- xref[1] -xref[2])
        end
    end
end 