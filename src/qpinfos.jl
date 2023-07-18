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