[![Build status](https://github.com/chmerdon/GradientRobustMultiPhysics.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/chmerdon/GradientRobustMultiPhysics.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://chmerdon.github.io/GradientRobustMultiPhysics.jl/stable/index.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://chmerdon.github.io/GradientRobustMultiPhysics.jl/dev/index.html)
[![DOI](https://zenodo.org/badge/229078096.svg)](https://zenodo.org/badge/latestdoi/229078096)

# GradientRobustMultiPhysics.jl

This package offers (mostly low-order) finite element methods for multiphysics problems in Julia that focus on the preservation of structural and qualitative properties, in particular the gradient-robustness property for the discretisation of (nearly) incompressible flows and resulting qualitative properties in coupled processes. The code therefore offers several classical and novel non-standard finite element discretisations to play and compare with in these applications and a toolkit to setup multi-physics problems by defining PDE systems and generating fixed-point iterations to solve them.

The implementation is based on [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl) that allows to have unstructured grids with mixed element geometries in it, e.g. triangles and quads in the same mesh.

Also note, that this package is part of the meta-package [PDELIB.jl](https://github.com/WIAS-BERLIN/PDELib.jl)

!!! note

    The focus is (at least currently) not on high-performance, high-order or parallel-computing. Also, this package is still in an early development stage and features and interfaces might change in future updates.
    

## Installation
via Julia package manager in Julia 1.6 or above:

```julia
# latest stable version
(@v1.6) pkg> add GradientRobustMultiPhysics
# latest version
(@v1.6) pkg> add GradientRobustMultiPhysics#master
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

Please have a look at the Examples.




## What is gradient-robustness?

Gradient-robustness is a feature of discretisations that exactly balance gradient forces in the momentum balance. In the case of the incompressible Navier--Stokes equations this means that the discrete velocity does not depend on the exact pressure. Divergence-free finite element methods have this property but are usually expensive and difficult to contruct. However, also non-divergence-free classical finite element methods can be made pressure-robust with the help of reconstruction operators applied to testfuntions in certain terms of the momentum balance, see e.g. references [1,2] below.

Recently gradient-robustness was also connected to the design of well-balanced schemes e.g. in the context of (nearly) compressible flows, see e.g. reference [3] below.

#### References

- [1]   "On the divergence constraint in mixed finite element methods for incompressible flows",\
        V. John, A. Linke, C. Merdon, M. Neilan and L. Rebholz,\
        SIAM Review 59(3) (2017), 492--544,\
        [>Journal-Link<](https://doi.org/10.1137/15M1047696),
        [>Preprint-Link<](http://www.wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=2015&number=2177)
- [2]   "Pressure-robustness and discrete Helmholtz projectors in mixed finite element methods for the incompressible Navier--Stokes equations",\
        A. Linke and C. Merdon,
        Computer Methods in Applied Mechanics and Engineering 311 (2016), 304--326,\
        [>Journal-Link<](http://dx.doi.org/10.1016/j.cma.2016.08.018)
        [>Preprint-Link<](http://www.wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=2016&number=2250)
- [3]   "A gradient-robust well-balanced scheme for the compressible isothermal Stokes problem",\
        M. Akbas, T. Gallouet, A. Gassmann, A. Linke and C. Merdon,\
        Computer Methods in Applied Mechanics and Engineering 367 (2020),\
        [>Journal-Link<](https://doi.org/10.1016/j.cma.2020.113069)
        [>Preprint-Link<](https://arxiv.org/abs/1911.01295)

