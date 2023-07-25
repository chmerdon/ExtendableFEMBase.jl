## object that shares information about the current quadrature point
mutable struct QPInfos{Ti, Tv, Ttime, Tx, Txref, TvG, TiG, PT}
    item::Ti
    cell::Ti
    region::Ti
    volume::Tv
    time::Ttime
    x::Vector{Tx}
    xref::Vector{Txref}
    grid::ExtendableGrid{TvG,TiG}
    params::PT
end

function QPInfos(xgrid::ExtendableGrid{Tv,Ti}; time = 0.0, dim = size(xgrid[Coordinates],1), T = Tv, x = ones(T, dim), params = []) where {Tv,Ti}
    return QPInfos(Ti(0), Ti(0), Ti(0), Tv(0), time, x, ones(T, dim), xgrid, params)
end


function standard_kernel(result, input, qpinfo)
    result .= input
    return nothing
end

function constant_one_kernel(result, qpinfo)
    result .= 1
    return nothing
end