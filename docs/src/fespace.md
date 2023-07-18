
# Finite Element Spaces and Arrays

This page describes the structure [FESpace](@ref) that acts as a finite element space on a given grid
and provides the associated degree of freedom maps [DofMaps](@ref) on demand.
See [Implemented Finite Elements](@ref) for a list of available finite element types.

Moreover, there are special arrays [FEVector](@ref) and [FEMatrix](@ref) that carry coefficients and discretised PDEOperators.


## FESpace

To generate a finite element space only a finite element type and a grid is needed, dofmaps are generated automatically on demand.

```@autodocs
Modules = [GradientRobustMultiPhysics]
Pages = ["finiteelements.jl"]
Order   = [:type, :function]
```

## DofMaps

```@autodocs
Modules = [GradientRobustMultiPhysics]
Pages = ["dofmaps.jl"]
Order   = [:type, :function]
```


The following DofMap subtypes are available and are used as keys to access the dofmap via ```FESpace[DofMap]``` (which is equivalent to ```FESpace.dofmaps[DofMap]```).

| DofMap             | Explanation                                       |
| :----------------: | :------------------------------------------------ | 
| CellDofs           | degrees of freedom for on each cell               | 
| FaceDofs           | degrees of freedom for each face                  | 
| EdgeDofs           | degrees of freedom for each edge (in 3D)          | 
| BFaceDofs          | degrees of freedom for each boundary face         |
| BEdgeDofs          | degrees of freedom for each boundary edge (in 3D) |


## FEVector

A FEVector consists of FEVectorBlocks that share a common one-dimensional array. Each block is associated to a FESpace and can only write into a region of the common array specified by offsets. It also acts as a one-dimensional AbstractArray itself.


```@autodocs
Modules = [GradientRobustMultiPhysics]
Pages = ["fevector.jl"]
Order   = [:type, :function]
```


## FEMatrix

A FEMatrix consists of FEMatrixBlocks that share a common ExtendableSparseMatrix. Each block is associated to two FESpaces and can only write into a submatrix of the common sparse matrix specified by offsets. It also acts as a two-dimensional AbstractArray itself.

```@autodocs
Modules = [GradientRobustMultiPhysics]
Pages = ["fematrix.jl"]
Order   = [:type, :function]
```