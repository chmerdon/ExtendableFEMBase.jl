################
### FEVector ###
################
#
# used to store coefficients for FESpaces and can have several blocks of different FESpaces
# acts like an AbstractArray{T,1}

"""
$(TYPEDEF)

block of an FEVector that carries coefficients for an associated FESpace and can be assigned as an AbstractArray (getindex, setindex, size, length)
"""
struct FEVectorBlock{T,Tv,Ti,FEType,APT} <: AbstractArray{T,1}
    name::String
    FES::FESpace{Tv,Ti,FEType,APT}
    offset::Int
    last_index::Int
    entries::Array{T,1} # shares with parent object
end


function Base.show(io::IO, FEB::FEVectorBlock; tol = 1e-14)
    @printf(io, "\n");
    #for j=1:FEB.offset+1:FEB.last_index
    #    if FEB.entries[j] > tol
    #        @printf(io, "[%d] = +%.3e", j, FEB.entries[j]);
    #    elseif FEB.entries[j] < -tol
    #        @printf(io, "[%d] %.3e", j, FEB.entries[j]);
    #    else
    #        @printf(io, "[%d] ********", j );
    #    end
    #    @printf(io, "\n");
    #end    
end

"""
$(TYPEDEF)

a plain array but with an additional layer of several FEVectorBlock subdivisions each carrying coefficients for their associated FESpace
"""
struct FEVector{T,Tv,Ti} #<: AbstractVector{T}
    FEVectorBlocks::Array{FEVectorBlock{T,Tv,Ti},1}
    entries::Array{T,1}
    tags::Vector{Any}
end

# overload stuff for AbstractArray{T,1} behaviour
Base.getindex(FEF::FEVector{T,Tv,Ti}, tag)  where {T,Tv,Ti} = FEF.FEVectorBlocks[findfirst(==(tag), FEF.tags)]
Base.getindex(FEF::FEVector,i::Int) = FEF.FEVectorBlocks[i]
Base.getindex(FEB::FEVectorBlock,i::Int)=FEB.entries[FEB.offset+i]
Base.getindex(FEB::FEVectorBlock,i::AbstractArray)=FEB.entries[FEB.offset.+i]
Base.getindex(FEB::FEVectorBlock,::Colon)=FEB.entries[FEB.offset+1:FEB.last_index]
Base.setindex!(FEB::FEVectorBlock, v, i::Int) = (FEB.entries[FEB.offset+i] = v)
Base.setindex!(FEB::FEVectorBlock, v, ::Colon) = (FEB.entries[FEB.offset+1:FEB.last_index] = v)
Base.setindex!(FEB::FEVectorBlock, v, i::AbstractArray) = (FEB.entries[FEB.offset.+i] = v)
Base.eltype(::FEVector{T}) where {T}  = T
Base.size(FEF::FEVector)=size(FEF.FEVectorBlocks)
Base.size(FEB::FEVectorBlock)=FEB.last_index-FEB.offset
Base.view(FEB::FEVectorBlock)=view(FEB.entries,FEB.offset+1:FEB.last_index)

function LinearAlgebra.norm(FEV::FEVector, p::Real = 2)
    return norm(FEV.entries)
end

function LinearAlgebra.norm(FEV::FEVectorBlock, p::Real = 2)
    return norm(view(FEV))
end

function norms(FEV::FEVector{T}, p::Real = 2) where {T}
    norms = zeros(T, length(FEV))
    for j = 1 : length(FEV)
        norms[j] = norm(view(FEV[j]), p)
    end
    return norms
end

get_ncomponents(FB::FEVectorBlock) = get_ncomponents(get_FEType(FB.FES))


"""
$(TYPEDSIGNATURES)

Returns the vector of FEspaces for the blocks of the given FEVector.
"""
function FESpaces(FEV::FEVector{T,Tv,Ti}) where {T,Tv,Ti}
    FEs = Array{FESpace{Tv,Ti},1}([])
    for j=1 : length(FEV.FEVectorBlocks)
        push!(FEs,FEV.FEVectorBlocks[j].FES)
    end   
    return FEs 
end


"""
$(TYPEDSIGNATURES)

Custom `length` function for `FEVector` that gives the number of defined FEMatrixBlocks in it
"""
Base.length(FEF::FEVector)=length(FEF.FEVectorBlocks)

"""
$(TYPEDSIGNATURES)

Custom `length` function for `FEVectorBlock` that gives the coressponding number of degrees of freedoms of the associated FESpace
"""
Base.length(FEB::FEVectorBlock)=FEB.last_index-FEB.offset

function FEVector(name::Union{String, Array{String,1}}, FES)
    return FEVector(FES; name = name)
end

"""
````
FEVector{T}(FES; name = "auto") where T <: Real
````

Creates FEVector that has one block if FES is a single FESpace, and a blockwise FEVector if FES is a vector of FESpaces.
Optionally a name for the vector (as a String), or each of the blocks (as a vector of Strings) can be specified.
"""
function FEVector(FES::FESpace{Tv,Ti,FEType,APT}; kwargs...) where {Tv,Ti,FEType,APT}
    return FEVector{Float64}([FES]; kwargs...)
end
function FEVector{T}(FES::FESpace{Tv,Ti,FEType,APT}; kwargs...) where {T,Tv,Ti,FEType,APT}
    return FEVector{T}([FES]; kwargs...)
end
function FEVector(FES::Array{<:FESpace{Tv,Ti},1}; kwargs...) where {Tv,Ti}
    return FEVector{Float64}(FES; kwargs...)
end
# main constructor
function FEVector{T}(FES::Array{<:FESpace{Tv,Ti},1}; name = "auto", tags = [], kwargs...) where {T,Tv,Ti}
    if name == "auto"
        names = ["#$j" for j in 1 : length(FES)]
    elseif typeof(name) == String
        names = Array{String,1}(undef, length(FES))
        for j = 1:length(FES)
            names[j] = name * " [#$j]"
        end    
    else
        names = name
    end
    @assert length(names) == length(FES)
    @debug "Creating FEVector mit blocks $((p->p.name).(FES))"
    ndofs = 0
    for j = 1:length(FES)
        ndofs += FES[j].ndofs
    end    
    entries = zeros(T,ndofs)
    Blocks = Array{FEVectorBlock{T,Tv,Ti},1}(undef,length(FES))
    offset = 0
    for j = 1:length(FES)
        Blocks[j] = FEVectorBlock{T,Tv,Ti,eltype(FES[j]),assemblytype(FES[j])}(names[j], FES[j], offset , offset+FES[j].ndofs, entries)
        offset += FES[j].ndofs
    end    
    return FEVector{T,Tv,Ti}(Blocks, entries, tags)
end


"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEVector` that prints some information on its blocks.
"""
function Base.show(io::IO, FEF::FEVector)
	println(io,"\nFEVector information")
    println(io,"====================")
    println(io,"   block  |  ndofs \t| FEType       \t\t| name/tag")
    for j=1:length(FEF)
        @printf(io," [%5d]  | ",j);
        @printf(io," %6d\t|",FEF[j].FES.ndofs);
        if length(FEF.tags) >= j
            @printf(io," %s  \t| %s/%s",FEF[j].FES.name,FEF[j].name,FEF.tags[j]);
        else
            @printf(io," %s  \t| %s",FEF[j].FES.name,FEF[j].name);
        end
    end    
end



"""
$(TYPEDSIGNATURES)

Custom `append` function for `FEVector` that adds a FEVectorBlock at the end.
"""
function Base.append!(FEF::FEVector{T}, FES::FESpace{Tv,Ti,FEType,APT}; name = "", tag = nothing) where {T,Tv,Ti,FEType,APT}
    append!(FEF.entries,zeros(T,FES.ndofs))
    newBlock = FEVectorBlock{T,Tv,Ti,FEType,APT}(name, FES, FEF.FEVectorBlocks[end].last_index , FEF.FEVectorBlocks[end].last_index+FES.ndofs, FEF.entries)
    push!(FEF.FEVectorBlocks,newBlock)
    if tag !== nothing
        push!(FEF.tags, tag)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Custom `fill` function for `FEVectorBlock` (only fills the block, not the complete FEVector).
"""
function Base.fill!(b::FEVectorBlock, value)
    for j = b.offset+1 : b.last_index
        b.entries[j] = value
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Adds FEVectorBlock b to FEVectorBlock a.
"""
function addblock!(a::FEVectorBlock, b::FEVectorBlock; factor = 1)
    l::Int = length(b)
    aoffset::Int = a.offset
    boffset::Int = b.offset
    for j = 1 : l
        a.entries[aoffset+j] += b.entries[boffset+j] * factor
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Adds Array b to FEVectorBlock a.
"""
function addblock!(a::FEVectorBlock, b::AbstractVector; factor = 1)
    aoffset::Int = a.offset
    for j = 1 : length(b)
        a.entries[aoffset+j] += b[j] * factor
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Scalar product between two FEVEctorBlocks
"""
function LinearAlgebra.dot(a::FEVectorBlock{T}, b::FEVectorBlock{T}) where {T}
    aoffset::Int = a.offset
    boffset::Int = b.offset
    result::T = 0.0
    for j = 1 : length(b)
       result += a.entries[aoffset+j] * b.entries[boffset+j]
    end
    return result
end
