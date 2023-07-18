# FEMatrix

A FEMatrix consists of FEMatrixBlocks that share a common ExtendableSparseMatrix. Each block is associated to two FESpaces and can only write into a submatrix of the common sparse matrix specified by offsets. It also acts as a two-dimensional AbstractArray itself.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["fematrix.jl"]
Order   = [:type, :function]
```