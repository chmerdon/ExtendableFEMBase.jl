
# Function Operators

## Available Operators

FunctionOperators are building blocks for the weak form and define the operations that should be applied to the trial and test functions (and their discrete representatives) inside some PDEOperator. Below is a list of currently available FunctionOperators. Note, that not all operators can be applied to all finite element types in principle, but potentially have to be understood in a broken sense and only make sense on certain parts of the mesh
(e.g. NormalFlux only on a face). Also note that all evaluations are returned as a vector,  so e.g.\ a gradient of a 2d vector-field will be a vector of length 4 (ordered component-wise).


| Function operator                                    | Description                                              | Mathematically                                                            |
| :--------------------------------------------------- | :------------------------------------------------------- | :------------------------------------------------------------------------ |
| Identity                                             | identity                                                 | ``v \rightarrow v``                                                       |
| IdentityComponent{c}                                 | identity of c-th component                               | ``v \rightarrow v_c``                                                     |
| NormalFlux                                           | normal flux (function times normal)                      | ``v \rightarrow v \cdot \vec{n}`` (only ON_FACES)                         |
| TangentFlux                                          | tangent flux (function times tangent)                    | ``v \rightarrow v \cdot \vec{t}`` (only ON_EDGES)                         |
| Gradient                                             | gradient/Jacobian (as a vector)                          | ``v \rightarrow \nabla v``                                                |
| SymmetricGradient                                    | symmetric part of the gradient                           | ``v \rightarrow Voigt(\mathrm{sym}(\nabla v))``                           |
| Divergence                                           | divergence                                               | ``v \rightarrow \mathrm{div}(v) = \nabla \cdot v``                        |
| CurlScalar                                           | curl operator 1D to 2D (rotated gradient)                | ``v \rightarrow [-dv/dx_2,dv/dx_1]``                                      |
| Curl2D                                               | curl operator 2D to 1D                                   | ``v \rightarrow dv_1/dx_2 - dv_2/dx_1``                                   |
| Curl3D                                               | curl operator 3D to 3D                                   | ``v \rightarrow \nabla \times v``                                         |
| Hessian                                              | Hesse matrix = all 2nd order derivatives (as a vector)   | ``v \rightarrow D^2 v``      (e.g. in 2D: xx,xy,yx,yy for each component) |
| SymmetricHessian{a}                                  | symmetric part of Hesse matrix, offdiagonals scaled by a | ``v \rightarrow sym(D^2 v)`` (e.g. in 2D: xx,yy,a*xy for each component)  |
| Laplacian                                            | Laplace Operator (diagonal of Hessian)                   | ``v \rightarrow \Delta v``   (e.g. in 2D: xx,yy for each component)       |



!!! note

    As each finite element type is transformed differently from the reference domain to the general domain,
    the evaluation of each function operator has to be implemented for each finite element class.
    Currently, not every function operator works in any dimension and for any finite element. More evaluations
    are added as soon as they are needed (and possibly upon request).
    Also, the function operators can be combined with user-defined actions to evaluate other operators that
    can be build from the ones available (e.g. the deviator).


## Jumps and Averages and Parents

If one of the operators above is evaluted ON_FACES for a finite element that is not continuous there, the code usually will crash or produce weird results. However, some operators can be transformed into a Jump- or Average operator and then either the jumps or the average of this operator along the face is assembled. The operator Jump(Identity) for example gives the jump of the identity evaluation on both sides of the face.
Seperate values of a discontinuous quantity on each neighbour of the face can be obtained using the Parent{1} and Parent{2} evaluations.

```@docs
Jump
Average
Parent
```

!!! note

    Currently this feature is only available for assembly on faces (2D and 3D) and certain function operators like Identity, Gradient, Laplacian, ReconstructionIdentity, ReconstructionGradient, NormalFlux, TangentFlux, but more
    are added as soon as they are needed (and possibly upon request).

    Also note that a Jump or Average operator has different behaviour depending on the Assembly Pattern it is used in. Usually, the input of the action used in the assembly pattern has the evaluation on one of the two neighbours at a time, which should be okay in a linear context. Only in ItemIntegrators the whole jump comes in the input (here the user can split the jump into left and right value via seperate Parent{1} and Parent{2} evaluations). In NonlinearForms jumps and averages should better not be used currently.


## Reconstruction Operators

There are special operators (see Table below) that allow to evaluate a usual operator of some discrete
reconstructed version of a vector-valued testfunction. These operators keep the discrete divergence exactly and so allow
for gradient-robust discretisations with classical non divergence-conforming ansatz spaces.
So far such operators are available for the vector-valued Crouzeix-Raviart and Bernardi--Raugel finite element types.


| Function operator                                    | Description                                             |
| :--------------------------------------------------- | :------------------------------------------------------ |
| ReconstructionIdentity{FEType}                       | reconstruction operator into specified FEType           |
| ReconstructionDivergence{FEType}                     | divergence of FEType reconstruction operator            |
| ReconstructionGradient{FEType}                       | gradient of FEType reconstruction operator              |


!!! note

    Currently this feature works with FEType = HdivRT0{d} and FEType = HdivBDM1{d} where d is the space dimension. However, solve! on a PDEDescription that includes these
    operators will only work if the function operators are at spots were it is applied to functions of type H1BR, H1CR or H1P2B.
    More reconstruction operators will be implemented at some later point.



## Operator Pairs (experimental)

Two function operators can be put into an OperatorPair so that one can provide effectively two operators in each argument of an assembly pattern. However, the user should make sure that both operators can be evaluated together reasonably (meaning both should be well-defined on the element geometries and the finite element space where the argument will be evaluated, and the action of the operator has to operate with coressponding input and result fields). This feature is still experimental and might have issues in some cases. OperatorTriple for a combination of three operators is also available.

```@docs
OperatorPair
OperatorTriple
```