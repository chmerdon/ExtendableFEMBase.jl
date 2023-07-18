
# FEEvaluator

FEEvaluators provide a structure that handles the evaluation of finite element basis functions for a given function operator, quadrature rule and item geometry. It stores the evaluations on the reference geometry (where derivatives are computed by automatic differentiation) and on the current mesh item. The current mesh item can be changed via the update! call.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["feevaluator.jl"]
Order   = [:type, :function]
```