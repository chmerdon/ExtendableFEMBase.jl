[![Build status](https://github.com/chmerdon/ExtendableFEMBase.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/chmerdon/ExtendableFEMBase.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://chmerdon.github.io/ExtendableFEMBase.jl/stable/index.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://chmerdon.github.io/ExtendableFEMBase.jl/dev/index.html)
[![DOI](https://zenodo.org/badge/229078096.svg)](https://zenodo.org/badge/latestdoi/229078096)

# ExtendableFEMBase.jl

This package offers (mostly low-order) finite element methods for multiphysics problems in Julia that focus on the preservation of structural and qualitative properties, in particular the gradient-robustness property for the discretisation of (nearly) incompressible flows and resulting qualitative properties in coupled processes. The code therefore offers several classical and novel non-standard finite element discretisations to play and compare with in these applications and a toolkit to setup multi-physics problems by defining PDE systems and generating fixed-point iterations to solve them.

The implementation is based on [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl) that allows to have unstructured grids with mixed element geometries in it, e.g. triangles and quads in the same mesh.

Also note, that this package is part of the meta-package [PDELIB.jl](https://github.com/WIAS-BERLIN/PDELib.jl)

!!! note

    The focus is (at least currently) not on high-performance, high-order or parallel-computing. Also, this package is still in an early development stage and features and interfaces might change in future updates.
    

## Installation
via Julia package manager in Julia 1.6 or above:

```julia
# latest stable version
(@v1.6) pkg> add ExtendableFEMBase
# latest version
(@v1.6) pkg> add ExtendableFEMBase#master
```

#### Dependencies on other Julia packages

[ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl)\
[GridVisualize.jl](https://github.com/j-fu/GridVisualize.jl)\
[ExtendableSparse.jl](https://github.com/j-fu/ExtendableSparse.jl)\
[DocStringExtensions.jl](https://github.com/JuliaDocs/DocStringExtensions.jl)\
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)\
[DiffResults.jl](https://github.com/JuliaDiff/DiffResults.jl)\



## Getting started

The general work-flow is as follows:

1. Mesh the domain of computation, possibly using one of the constructors by ExtendableGrid.jl or via mesh generators in [SimplexGridFactory.jl](https://github.com/j-fu/SimplexGridFactory.jl).
2. Describe your PDE system with the help of the [PDE Description](@ref) and [PDE Operators](@ref). User parameters and customised operator actions are framed with the help of [User Data and Actions](@ref).
3. Discretise, i.e. choose suitable finite element ansatz spaces for the unknowns of your PDE system.
4. Solve (stationary, time-dependent, iteratively?)
5. Postprocess (compute stuff, plot, export data)

Please have a look at the Pluto notebook examples.