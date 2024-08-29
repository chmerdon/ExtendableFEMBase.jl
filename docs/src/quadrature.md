
# Quadrature

Usually quadrature is a hidden layer as quadrature rules are chosen automatically based on the polynomial degree of the ansatz functions and the specified quadorder of the user data.

Hence, quadrature rules are only needed if the user wants write his own low-level assembly.


Quadrature rules consist of points (coordinates of evaluation points with respect to reference geometry) and weights. There are constructors for several AbstractElementGeometries (from ExtendableGrids) and different order (some have generic formulas for arbitrary order), see below for a detailed list.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["quadrature.jl"]
Order   = [:type, :function]
```


#### Accumulating Vector (not relevant for users, but for completeness)

Internally a global integration uses an accumulating vector and calls the cell-wise integration.

```@docs
AccumulatingVector
```