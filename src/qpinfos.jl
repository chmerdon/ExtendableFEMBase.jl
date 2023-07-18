## object that shares information about the current quadrature point
mutable struct QPInfos{Ti, Tv, Ttime, Tx, Txref, TvG, TiG, PT}
    item::Ti
    region::Ti
    volume::Tv
    time::Ttime
    x::Vector{Tx}
    xref::Vector{Txref}
    grid::ExtendableGrid{TvG,TiG}
    params::PT
end

function QPInfos(xgrid::ExtendableGrid{Tv,Ti}; time = 0.0, dim = size(xgrid[Coordinates],1), T = Tv, params = []) where {Tv,Ti}
    return QPInfos(Ti(0), Ti(0), Tv(0), time, zeros(T, dim), zeros(T, dim), xgrid, params)
end