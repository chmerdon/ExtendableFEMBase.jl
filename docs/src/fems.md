
# Implemented Finite Elements

This page describes the finite element type-tree and lists all implemented finite elements.



## The Finite Element Type-Tree

Finite elements are abstract type leaves in a type-tree. The complete tree looks like this:

```
AbstractFiniteElement
├─ AbstractH1FiniteElement
│  ├─ AbstractH1FiniteElementWithCoefficients
│  │  ├─ H1P1TEB
│  │  └─ H1BR
│  ├─ H1CR
│  ├─ H1MINI
│  ├─ L2P0
│  ├─ L2P1
│  ├─ H1P1
│  ├─ H1P2
│  ├─ H1P2B
│  ├─ H1P3
│  ├─ H1Pk
│  ├─ H1Q1
│  └─ H1Q2
├─ AbstractHcurlFiniteElement
│  └─ HCURLN0
└─ AbstractHdivFiniteElement
   ├─ HDIVBDM1
   ├─ HDIVBDM2
   ├─ HDIVRT0
   └─ HDIVRT1
```


#### Remarks
- each type depends on one/two or three parameters, the first one is always the number of components (ncomponents) that determines if the
  finite element is scalar- or veector-valued; some elements additionaly require the parameter edim <: Int if they are structurally different in different space dimensions; arbitrary order elements require a third parameter that determines the order
- each finite elements mainly comes with a set of basis functions in reference coordinates for each applicable AbstractElementGeometry and degrees of freedom maps for each [Assembly Type](@ref) (coded as a string)
- broken finite elements are possible via the broken switch in the [FESpace](@ref) constructor
- the type steers how the basis functions are transformed from local to global coordinates and how FunctionOperators are evaluated
- depending on additional continuity properties of the element types more basis function sets are defined:
    - AbstractH1FiniteElements additionally have evaluations of nonzero basisfunctions on faces/bfaces
    - AbstractHdivFiniteElements additionally have evaluations of nonzero normalfluxes of basisfunctions on faces/bfaces
    - AbstractHcurlFiniteElements additionally have evaluations of nonzero tangentfluxes of basisfunctions on edges/bedges
- each finite element has its own implemented standard interpolation interpolate! (see [Finite Element Interpolations](@ref)) that can be applied to a [Data Function](@ref), below it is shortly described what this means for each finite element


## List of implemented Finite Elements

The following table lists all curently implemented finite elements and on which geometries they are available (in brackets a dofmap pattern for CellDofs is shown and the number of local degrees of freedom for a vector-valued realisation). Click on the FEType to find out more details.

| FEType | Triangle2D | Parallelogram2D | Tetrahedron3D | Parallelepiped3D |
| :----------------: | :----------------: |  :----------------: |  :----------------: |  :----------------: | 
| AbstractH1FiniteElementWithCoefficients |   |   |   |   |
| [`H1BR`](@ref) | ✓ (N1f1, 9) | ✓ (N1f1, 12) | ✓ (N1f1, 16) |   |
| [`H1P1TEB`](@ref) | ✓ (N1f1, 9) |   | ✓ (N1e1, 18) |   |
| AbstractH1FiniteElement |   |   |   |   |
| [`H1BUBBLE`](@ref) | ✓ (I1, 2) | ✓ (I1, 2) | ✓ (I1, 3) |   |
| [`H1CR`](@ref) | ✓ (F1, 6) | ✓ (F1, 8) | ✓ (F1, 12) |   |
| [`H1MINI`](@ref) | ✓ (N1I1, 8) | ✓ (N1I1, 10) | ✓ (N1I1, 15) |   |
| [`L2P0`](@ref) | ✓ (I1, 2) | ✓ (I1, 2) | ✓ (I1, 3) | ✓ (I1, 3) |
| [`L2P1`](@ref) | ✓ (I3, 6) | ✓ (I3, 6) | ✓ (I4, 12) | ✓ (I4, 12) |
| [`H1P1`](@ref) | ✓ (N1, 6) |  | ✓ (N1, 12) |  |
| [`H1P2`](@ref) | ✓ (N1F1, 12) |  | ✓ (N1E1, 30) |   |
| [`H1P2B`](@ref) | ✓ (N1F1I1, 14) |   |   |   |
| [`H1P3`](@ref) | ✓ (N1F2I1, 20) |   | ✓ (N1E2F1, 60)  |   |
| [`H1Pk`](@ref) | ✓ (order-dep) |   |   |   |
| [`H1Q1`](@ref) | ✓ (N1, 6) | ✓ (N1, 8) | ✓ (N1, 12) | ✓ (N1, 24) |
| [`H1Q2`](@ref) | ✓ (N1F1, 12) | ✓ (N1F1I1, 18) | ✓ (N1E1, 30) |   |
| AbstractHcurlFiniteElement |   |   |   |   |
| [`HCURLN0`](@ref) | ✓ (f1, 3) | ✓ (f1, 4) | ✓ (e1, 6) |   |
| AbstractHdivFiniteElement |   |   |   |   |
| [`HDIVBDM1`](@ref) | ✓ (f2, 6) | ✓ (f2, 8) | ✓ (f3, 12) |   |
| [`HDIVBDM2`](@ref) | ✓ (f3i3, 12) |   |   |   |
| [`HDIVRT0`](@ref) | ✓ (f1, 3) | ✓ (f1, 4) | ✓ (f1, 4) | ✓ (f1, 6) |
| [`HDIVRT1`](@ref) | ✓ (f2i2, 8) |   | ✓ (f3i3, 15) |   |


Note: the dofmap pattern describes the connection of the local degrees of freedom to entities of the grid and also hints to the continuity. Here, "N" or "n" means nodes, "F" or "f" means faces, "E" or "e" means edges and "I" means interior (dofs without any continuity across elements). Capital letters cause that every component has its own degree of freedom, while small letters signalize that only one dof is associated to the entity. As an example "N1f1" (for the Bernardi-Raugel element) means that at each node sits one dof per component and at each face sits a single dof. Usually finite elements that involve small letters are only defined vector-valued (i.e. the number of components has to match the element dimension), while finite elements that only involve capital letters are available for any number of components.


## H1-conforming finite elements

### P0 finite element

Piecewise constant finite element that has one degree of freedom on each cell of the grid. (It is masked as a H1-conforming finite element, because it uses the same operator evaulations.)

The interpolation of a given function into this space preserves the cell integrals.

```@docs
L2P0
```

### P1 finite element

The lowest-order Courant finite element that has a degree of freedom on each vertex of the grid. On simplices the
basis functions coincide with the linear barycentric coordinates. Only the L2P1 element is also defined on quads.

The interpolation of a given function into this space performs point evaluations at the nodes.

```@docs
L2P1
H1P1
```

### Q1 finite element

The lowest-order finite element that has a degree of freedom on each vertex of the grid. On simplices the
basis functions coincide with the linear barycentric coordinates. This element is also defined on quads.

The interpolation of a given function into this space performs point evaluations at the nodes.

```@docs
H1Q1
```


### MINI finite element

The mini finite element adds cell bubles to the P1 element that are e.g. beneficial to define inf-sup stable finite element pairs for the Stokes problem.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves its
cell integral.

```@docs
H1MINI
```


### P1TEB finite element

This element adds tangent-weighted edge bubbles to the P1 finite element and therefore is only available as a vector-valued element.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves face integrals of its tangential flux.

```@docs
H1P1TEB
```


### Bernardi-Raugel (BR) finite element

The Bernardi-Raugel adds normal-weighted face bubbles to the P1 finite element and therefore is only available
as a vector-valued element.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves face integrals of its normal flux.

```@docs
H1BR
```


### P2 finite element

The P2 finite element method on simplices equals quadratic polynomials. On the Triangle2D shape the degrees of freedom
are associated with the three vertices and the three faces of the triangle. On the Tetrahedron3D shape the degrees of freedom are associated with the four verties and the six edges.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves its face/edge integrals in 2D/3D.

```@docs
H1P2
```


### Q2 finite element

A second order finite element. On simplices it equals the P2 finite element, and on Quadrilateral2D it has 9 degrees of freedom (vertices, faces and one cell bubble).

The interpolation of a given function into this space performs point evaluations at the nodes and preserves lowest order face moments and (only on quads) also the cell integreal mean.

```@docs
H1Q2
```

### P2B finite element

The P2B finite element adds additional cell bubles (in 2D and 3D) and face bubbles (only in 3D) that are e.g. used to define inf-sup stable finite element pairs for the Stokes problem.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves its cell and face integrals in 2D and also edge integrals in 3D.

```@docs
H1P2B
```

### P3 finite element

The P3 finite element method on simplices equals cubic polynomials. On the Triangle2D shape the degrees of freedom
are associated with the three vertices, the three faces (double dof) of the triangle and the cell itself (one cell bubble).

The interpolation of a given function into this space performs point evaluations at the nodes and preserves cell and face integrals in 2D.

```@docs
H1P3
```

### Pk finite element (experimental)

The Pk finite element method generically generates polynomials of abitrary order k on simplices (Edge1D, Triangle2D so far).

The interpolation of a given function into this space performs point evaluations at the nodes and preserves cell and face integrals in 2D (moment order depends on the order and the element dimension).

```@docs
H1Pk
```

### Crouzeix-Raviart (CR) finite element

The Crouzeix-Raviart element associates one lowest-order function with each face. On the Triangle2D shape, the basis function of a face is one minus two times the nodal basis function of the opposite node. 

The interpolation of a given function into this space preserves its face integrals.

```@docs
H1CR
```



## Hdiv-conforming finite elements

These Raviart-Thomas and Brezzi-Douglas-Marini finite elements of lower order and their standard interpolations are available:

```@docs
HDIVRT0
HDIVBDM1
HDIVRT1
HDIVBDM2
```

## Hcurl-conforming finite elements

So far only the lowest order Nedelec element is available in 2D and 3D. On Triangle2D it has one degree of freedom for each face (i.e. the rotated RT0 element), on Tetrahedron3D it has one degree of freedom associated to each of the six edges.

Its standard interpolation of a given functions preserves its tangential face/edge integrals.

```@docs
HCURLN0
```

