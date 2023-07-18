"""
````
abstract type HDIVRT1BUBBLE{edim} <: AbstractHdivFiniteElement where {edim<:Int}
````

Just internal (normal-zero) Hdiv-conforming vector-valued (ncomponents = edim) Raviart-Thomas space of order 1.

allowed ElementGeometries:
- Triangle2D
"""
abstract type HDIVRT1BUBBLE{edim} <: AbstractHdivFiniteElement where {edim<:Int} end

function Base.show(io::Core.IO, FEType::Type{<:HDIVRT1BUBBLE{edim}}) where {edim}
    print(io,"HDIVRT1BUBBLE{$edim}")
end

get_ncomponents(FEType::Type{<:HDIVRT1BUBBLE}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRT1BUBBLE}, EG::Type{<:AbstractElementGeometry1D}) = 0
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRT1BUBBLE}, EG::Type{<:Triangle2D}) = 0
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRT1BUBBLE}, EG::Type{<:Triangle2D}) = 2
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRT1BUBBLE}, EG::Type{<:Tetrahedron3D}) = 3

get_polynomialorder(::Type{<:HDIVRT1BUBBLE{2}}, ::Type{<:AbstractElementGeometry1D}) = 1;
get_polynomialorder(::Type{<:HDIVRT1BUBBLE{3}}, ::Type{<:AbstractElementGeometry2D}) = 1;
get_polynomialorder(::Type{<:HDIVRT1BUBBLE{2}}, ::Type{<:AbstractElementGeometry2D}) = 2;
get_polynomialorder(::Type{<:HDIVRT1BUBBLE{3}}, ::Type{<:AbstractElementGeometry3D}) = 2;

get_dofmap_pattern(FEType::Type{<:HDIVRT1BUBBLE{2}}, ::Type{CellDofs}, EG::Type{<:Triangle2D}) = "i2"
get_dofmap_pattern(FEType::Type{<:HDIVRT1BUBBLE{3}}, ::Type{CellDofs}, EG::Type{<:Tetrahedron3D}) = "i3"

isdefined(FEType::Type{<:HDIVRT1BUBBLE}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:HDIVRT1BUBBLE}, ::Type{<:Tetrahedron3D}) = true

function ExtendableGrids.interpolate!(Target::AbstractArray{T,1}, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_FACES}, exact_function!; items = [], time = 0) where {T,Tv,Ti,FEType <: HDIVRT1BUBBLE,APT}
    ncomponents = get_ncomponents(FEType)
    if items == []
        items = 1 : num_sources(FE.xgrid[FaceNodes])
    end

   # integrate normal flux of exact_function over edges
   xFaceNormals::Array{Tv,2} = FE.xgrid[FaceNormals]
   nfaces = num_sources(xFaceNormals)
   function normalflux_eval()
       temp = zeros(T,ncomponents)
       function closure(result, x, face)
            eval_data!(temp, exact_function!, x, time)
            result[1] = 0
            for j = 1 : ncomponents
               result[1] += temp[j] * xFaceNormals[j,face]
            end 
       end   
   end   
   edata_function = ExtendedDataFunction(normalflux_eval(), [1, ncomponents]; dependencies = "XI", bonus_quadorder = exact_function!.quadorder)
   integrate!(Target, FE.xgrid, ON_FACES, edata_function; items = items)
   
   # integrate first moment of normal flux of exact_function over edges
   function normalflux2_eval()
       temp = zeros(T,ncomponents)
       function closure(result, x, face, xref)
            eval_data!(temp, exact_function!, x, time)
            result[1] = 0.0
            for j = 1 : ncomponents
               result[1] += temp[j] * xFaceNormals[j,face]
            end
            result[1] *= (xref[1] - 1//ncomponents)
       end   
   end   
   edata_function2 = ExtendedDataFunction(normalflux2_eval(), [1, ncomponents]; dependencies = "XIL", bonus_quadorder = exact_function!.quadorder + 1)
   integrate!(Target, FE.xgrid, ON_FACES, edata_function2; items = items, index_offset = nfaces)

    if ncomponents == 3
        function normalflux3_eval()
            temp = zeros(T,ncomponents)
            function closure(result, x, face, xref)
                eval_data!(temp, exact_function!, x, time)
                result[1] = 0.0
                for j = 1 : ncomponents
                    result[1] += temp[j] * xFaceNormals[j,face]
                end
                result[1] *= (xref[2] - 1//ncomponents)
            end   
        end   
        edata_function3 = ExtendedDataFunction(normalflux3_eval(), [1, ncomponents]; dependencies = "XIL", bonus_quadorder = exact_function!.quadorder + 1)
        integrate!(Target, FE.xgrid, ON_FACES, edata_function3; items = items, time = time, index_offset = 2*nfaces)
    end
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T,1}, FE::FESpace{Tv,Ti,FEType,APT}, ::Type{ON_CELLS}, exact_function!; items = [], time = 0) where {T,Tv,Ti,FEType <: HDIVRT1BUBBLE,APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.xgrid[CellFaces], items)
    interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems)

    # set values of interior RT1 functions by integrating over cell
    # they are chosen such that integral mean of exact function is preserved on each cell
    ncomponents = get_ncomponents(FEType)
    ncells = num_sources(FE.xgrid[CellNodes])
    xCellVolumes::Array{Tv,1} = FE.xgrid[CellVolumes]
    xCellDofs::DofMapTypes{Ti} = FE[CellDofs]
    means = zeros(T,ncomponents,ncells)
    integrate!(means, FE.xgrid, ON_CELLS, exact_function!)
    EG = (ncomponents == 2) ? Triangle2D : Tetrahedron3D
    qf = QuadratureRule{T,EG}(2)
    FEB = FEBasisEvaluator{T,EG,Identity,ON_CELLS}(FE, qf)
    if items == []
        items = 1 : ncells
    end

    basisval = zeros(T,ncomponents)
    IMM = zeros(T,ncomponents,ncomponents)
    interiordofs = zeros(Int,ncomponents)
    interior_offset::Int = (ncomponents == 2) ? 6 : 12
    for cell in items
        update_basis!(FEB,cell)
        # compute mean value of facial RT1 dofs
        for dof = 1 : interior_offset
            for i = 1 : length(qf.w)
                eval_febe!(basisval,FEB, dof, i)
                for k = 1 : ncomponents
                    means[k,cell] -= basisval[k] * Target[xCellDofs[dof,cell]] * xCellVolumes[cell] * qf.w[i]
                end
            end
        end
        # compute mss matrix of interior dofs
        fill!(IMM,0)
        for dof = 1:ncomponents
            for i = 1 : length(qf.w)
                eval_febe!(basisval,FEB, interior_offset + dof, i)
                for k = 1 : ncomponents
                    IMM[k,dof] += basisval[k] * xCellVolumes[cell] * qf.w[i]
                end
            end
            interiordofs[dof] = xCellDofs[interior_offset + dof,cell] 
        end
        Target[interiordofs] = IMM\means[:,cell]
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT1BUBBLE{2}}, ::Type{<:Triangle2D})
    function closure(refbasis,xref)
        # RT0 basis
        refbasis[1,1] = 12*xref[2] * xref[1];        refbasis[1,2] = 12*xref[2] * (xref[2]-1)
        refbasis[2,1] = 12*xref[1] * (xref[1]-1);    refbasis[2,2] = 12*xref[1] * xref[2]
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT1BUBBLE{3}}, ::Type{<:Tetrahedron3D})
    function closure(refbasis,xref)
        refbasis[end] = 1 - xref[1] - xref[2] - xref[3]
        # RT0 basis
        refbasis[1,1] = 2*xref[1];      refbasis[1,2] = 2*xref[2];      refbasis[1,3] = 2*(xref[3]-1)
        refbasis[2,1] = 2*xref[1];      refbasis[2,2] = 2*(xref[2]-1);  refbasis[2,3] = 2*xref[3]
        refbasis[3,1] = 2*xref[1];      refbasis[3,2] = 2*xref[2];      refbasis[3,3] = 2*xref[3]

        for k = 1 : 3
            # interior functions
            refbasis[1,k] = 12*xref[3] * refbasis[1,k];
            refbasis[2,k] = 12*xref[2] * refbasis[2,k];
            refbasis[3,k] = 12*xref[1] * refbasis[3,k];
        end
    end
end
