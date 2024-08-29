###########################
### ACCUMULATING VECTOR ###
###########################

"""
$(TYPEDEF)

vector that is acting as an AbstractArray{T, 2} and 
automatically accumulates all values from the second dimension

AV[k,j] += s for any j results in AV.entries[k] += s
"""
struct AccumulatingVector{T} <: AbstractArray{T, 2}
	entries::Array{T, 1}
	size2::Int
end

# 

# overload stuff for AbstractArray{T,2} behaviour
Base.getindex(AV::AccumulatingVector, i::Int, j) = AV.entries[i]
Base.getindex(AV::AccumulatingVector, i::AbstractArray, j) = AV.entries[i]
Base.getindex(AV::AccumulatingVector, ::Colon, j) = AV.entries
Base.setindex!(AV::AccumulatingVector, v, i::Int, j) = (AV.entries[i] = v)
Base.setindex!(AV::AccumulatingVector, v, ::Colon, j) = (AV.entries .= v)
Base.setindex!(AV::AccumulatingVector, v, i::AbstractArray, j) = (AV.entries[i] .= v)
Base.size(AV::AccumulatingVector) = [length(AV.entries), AV.size2]
Base.length(AV::AccumulatingVector) = length(AV.entries)

Base.fill!(AV::AccumulatingVector, v) = (fill!(AV.entries, v))
