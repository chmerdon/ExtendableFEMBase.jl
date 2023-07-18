###########################
### ACCUMULATING VECTOR ###
###########################
#
# used to store coefficients for FESpaces and can have several blocks of different FESpaces
# acts like an AbstractArray{T,1}

"""
$(TYPEDEF)

block of an FEVector that carries coefficients for an associated FESpace and can be assigned as an AbstractArray (getindex, setindex, size, length)
"""
struct AccumulatingVector{T} <: AbstractArray{T,2}
    entries::Array{T,1}
    size2::Int
end

# AV[k,j] += s for any j results in entries[k] += s

# overload stuff for AbstractArray{T,2} behaviour
Base.getindex(AV::AccumulatingVector,i::Int, j)=AV.entries[i]
Base.getindex(AV::AccumulatingVector,i::AbstractArray, j)=AV.entries[i]
Base.getindex(AV::AccumulatingVector,::Colon, j)=AV.entries
Base.setindex!(AV::AccumulatingVector, v, i::Int, j) = (AV.entries[i] = v)
Base.setindex!(AV::AccumulatingVector, v, ::Colon, j) = (AV.entries .= v)
Base.setindex!(AV::AccumulatingVector, v, i::AbstractArray, j) = (AV.entries[i] .= v)
Base.size(AV::AccumulatingVector)=[length(AV.entries),AV.size2]
Base.length(AV::AccumulatingVector)=length(AV.entries)

Base.fill!(AV::AccumulatingVector, v) = (fill!(AV.entries, v))