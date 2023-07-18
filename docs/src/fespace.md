
# FESpace

To generate a finite element space only a finite element type and a grid is needed, dofmaps are generated automatically on their first demand.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["finiteelements.jl"]
Order   = [:type, :function]
```

## DofMaps

```@autodocs
Modules = [ExtendableFEMBase]
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

