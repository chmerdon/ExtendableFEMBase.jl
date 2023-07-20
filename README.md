# ExtendableFEMBase

This package provides basic finite element structures to setup finite element schemes on ExtendableGrids:

- FEEvaluator (finite element basis evaluators for several H1, Hdiv and Hcurl elements)
- FESpace (Discrete finite element space with respect to a mesh from ExtendableGrids, knows the Dofmaps)
- FEMatrix (block overlay for an ExtendableSparse matrix, where each block corresponds to a coupling between two FESpaces in a system)
- FEVector (block overlay for an array, where each block corresponds to a FESpace)
- QuadratureRule (basic quadrature rules for different ElementGeometries from ExtendableGrids)
- interpolations (standard interpolations into the provided finite element spaces, averaging routines and interpolations between meshes/FESpaces)
- FunctionOperators (primitive linear operators like Identity, Gradient, Divergence) and rules how to evaluate them for for different finite element types
- reconstruction operators (special FunctionOperators that involve an interpolation into a different finite element type)
  
  
