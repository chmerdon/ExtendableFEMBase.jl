
# Finite Element Interpolations

## Source functions and QPInfo

The functions that can be interpolated with the methods below are expected to have a certain interface, i.e.:
```julia
function f!(result, qpinfo) end
```
The qpinfo argument communicates vast information of the current quadrature point:

| qpinfo child       | Type               | Description         |
| :----------------  | :----------------  |  :---------------- |
| qpinfo.x           | Vector{Real}       | space coordinates of quadrature point |
| qpinfo.time        | Real               | current time |
| qpinfo.item        | Integer            | current item that contains qpinfo.x |
| qpinfo.region      | Integer            | region number of item |
| qpinfo.xref        | Vector{Real}       | reference coordinates within item of qpinfo.x |
| qpinfo.volume      | Real               | volume of item |
| qpinfo.params      | Vector{Any}        | parameters that can be transfered via keyword arguments |


## Standard Interpolations

Each finite element has its standard interpolator that can be applied to some user-defined DataFunction. Instead of interpolating on the full cells, the interpolation can be restricted to faces or edges. 

It is also possible to interpolate finite element functions on one grid onto a finite element function on another grid (experimental feature, does not work for all finite elements yet and shall be extended to interpolations of operator evaluations as well in future).

```@docs
interpolate!
```

## Nodal Evaluations

Usually, Plotters need nodal values, so there is a gengeric function that evaluates any finite element function at the nodes of the grids (possibly by averaging if discontinuous). In case of Identity evaluations of an H1-conforming finite element, the function nodevalues_view can generate a view into the coefficient field that avoids further allocations.


```@docs
nodevalues!
nodevalues
nodevalues_view
```

## Lazy Interpolation

To interpolate between different finite element spaces and meshes, there is a lazy interpolation routine that
works in all cases (but is not very efficient as it involves a PointeEvaluator and CellFinder):

```@docs
lazy_interpolate!
```