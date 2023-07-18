# FEVector

A FEVector consists of FEVectorBlocks that share a common one-dimensional array. Each block is associated to a FESpace and can only write into a region of the common array specified by offsets. It also acts as a one-dimensional AbstractArray itself.


```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["fevector.jl"]
Order   = [:type, :function]
```
