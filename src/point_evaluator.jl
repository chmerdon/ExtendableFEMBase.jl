mutable struct PointEvaluator{Tv <: Real, UT, KFT <: Function}
	u_args::Array{UT, 1}
	ops_args::Array{DataType, 1}
	kernel::KFT
	BE_args::Any
	L2G::Any
	CF::Any
	lastitem::Int
	eval_selector::Any
	evaluator_bary::Any
	evaluator::Any
	xref::Vector{Tv}
	parameters::Dict{Symbol, Any}
end

default_peval_kwargs() = Dict{Symbol, Tuple{Any, String}}(
	:name => ("PointEvaluator", "name for operator used in printouts"),
	:resultdim => (0, "dimension of result field (default = length of operators)"),
	:params => (nothing, "array of parameters that should be made available in qpinfo argument of kernel function"),
	:verbosity => (0, "verbosity level"),
)


"""
````
function Pointevaluator(
	[kernel!::Function],
	oa_args::Array{<:Tuple{<:Any, DataType},1};
	kwargs...)
````

Generates a PointEvaluator that can evaluate the specified operator evaluations
at arbitrary points. If no kernel function is given, the arguments
are given directly. If a kernel is provided, the arguments are postprocessed
accordingly and the kernel has to be conform to the interface

	kernel!(result, eval_args, qpinfo)

where qpinfo allows to access information at the current evaluation point.
Additionally the length of the result needs to be specified via the kwargs.

Evaluation can be triggered via the evaluate function after an initialize! call.

Operator evaluations are tuples that pair a tag (to identify an unknown or the position in the vector)
with a FunctionOperator.

Keyword arguments:
$(_myprint(default_peval_kwargs()))
"""
function PointEvaluator(kernel, u_args, ops_args, sol = nothing; Tv = Float64, kwargs...)
	parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_peval_kwargs())
	_update_params!(parameters, kwargs)
	@assert length(u_args) == length(ops_args)
	PE = PointEvaluator{Tv, typeof(u_args[1]), typeof(kernel)}(u_args, ops_args, kernel, nothing, nothing, nothing, 1, nothing, nothing, nothing, zeros(Tv, 2), parameters)
	if sol !== nothing
		initialize!(PE, sol)
	end
	return PE
end

function PointEvaluator(kernel::Function, oa_args::Array{<:Tuple{<:Any, DataType}, 1}, sol = nothing; kwargs...)
	u_args = [oa[1] for oa in oa_args]
	ops_args = [oa[2] for oa in oa_args]
	return PointEvaluator(kernel, u_args, ops_args, sol; kwargs...)
end

function PointEvaluator(oa_args::Array{<:Tuple{<:Any, DataType}, 1}, sol = nothing; kwargs...)
	return PointEvaluator(standard_kernel, oa_args, sol; kwargs...)
end


"""
````
function initialize!(
	O::PointEvaluator,
	sol;
	time = 0,
	kwargs...)
````

Initializes the PointEvaluator for the specified solution.
"""
function initialize!(O::PointEvaluator{T, UT}, sol; time = 0, kwargs...) where {T, UT}
	_update_params!(O.parameters, kwargs)
	if UT <: Integer
		ind_args = O.u_args
	else
		ind_args = [findfirst(==(u), sol.tags) for u in O.u_args]
	end
	FES_args = [sol[j].FES for j in ind_args]
	nargs = length(FES_args)
	xgrid = FES_args[1].xgrid
	Ti = eltype(xgrid[CellNodes])
	EGs = xgrid[UniqueCellGeometries]
	AT = ON_CELLS
	gridAT = ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_args[1]), AT)
	xgrid = FES_args[1].xgrid
	itemregions = xgrid[CellRegions]
	itemgeometries = xgrid[CellGeometries]

	O.CF = CellFinder(xgrid)
	O.xref = zeros(T, size(xgrid[Coordinates],1))

	O.BE_args = Array{Array{<:FEEvaluator, 1}, 1}([])
	O.L2G = []
	for EG in EGs
		## FE basis evaluator for EG
		push!(O.BE_args, [FEEvaluator(FES_args[j], O.ops_args[j], QuadratureRule{T, EG}(0); AT = AT) for j in 1:nargs])

		## L2G map for EG
		push!(O.L2G, L2GTransformer(EG, xgrid, gridAT))
	end

	## parameter structure
	QPinfo = QPInfos(xgrid; time = time, params = O.parameters[:params])

	## prepare input args
	op_lengths_args = [size(O.BE_args[1][j].cvals, 1) for j ∈ 1:nargs]
	op_offsets_args = [0]
	append!(op_offsets_args, cumsum(op_lengths_args))
	input_args = zeros(T, op_offsets_args[end])

	FEATs_args = [ExtendableFEMBase.EffAT4AssemblyType(get_AT(FES_args[j]), AT) for j ∈ 1:nargs]
	itemdofs_args::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [FES_args[j][Dofmap4AssemblyType(FEATs_args[j])] for j ∈ 1:nargs]
	kernel = O.kernel

	function eval_selector(item)
		return findfirst(==(itemgeometries[item]), EGs)
	end

	function _evaluate_bary!(
		result,
		BE_args::Array{<:FEEvaluator, 1},
		L2G::L2GTransformer,
		xref,
		item, # cell used to evaluate local coordinates
	)

		for id ∈ 1:nargs
			# update basis evaluations at xref
			ExtendableFEMBase.relocate_xref!(BE_args[id], xref)

			# update operator eveluation on item
			update_basis!(BE_args[id], item)
		end

		# update QPinfo
		QPinfo.item = item
		QPinfo.region = itemregions[item]
		QPinfo.xref = xref
		update_trafo!(L2G, item)
		eval_trafo!(QPinfo.x, L2G, xref)
		# evaluate operator
		fill!(input_args, 0)
		for id ∈ 1:nargs
			for j ∈ 1:size(BE_args[id].cvals, 2)
				dof_j = itemdofs_args[id][j, item]
				for d ∈ 1:op_lengths_args[id]
					input_args[d+op_offsets_args[id]] += sol[ind_args[id]][dof_j] * BE_args[id].cvals[d, j, 1]
				end
			end
		end

		## evaluate kernel
		kernel(result, input_args, QPinfo)

		return nothing
	end

	## initialize cell finder
	CF = CellFinder(xgrid)
	xref = zeros(T, size(xgrid[Coordinates],1))
	function _evaluate!(
		result,
		BE_args::Array{<:FEEvaluator, 1},
		L2G::L2GTransformer,
		x
	)


		# evaluate in barycentric coordinates
		_evaluate_bary!(result, BE_args, L2G, xref, item)

		return nothing
	end
	O.evaluator = _evaluate!
	O.evaluator_bary = _evaluate_bary!
	O.eval_selector = eval_selector

	return nothing
end


"""
````
function evaluate_bary!(
	result,
	PE::PointEvaluator,
	xref, 
	item
	)
````

Evaluates the PointEvaluator at the specified reference coordinates in the cell with the specified item number.
"""
function evaluate_bary!(
	result,
	PE::PointEvaluator,
	xref,
	item,
)

	## find cell geometry id
	j = PE.eval_selector(item)

	## evaluate
	PE.evaluator_bary(result, PE.BE_args[j], PE.L2G[j], xref, item)
end

"""
````
function evaluate!(
	result,
	PE::PointEvaluator,
	x
	)
````

Evaluates the PointEvaluator at the specified coordinates x.
(To do so it internally calls CellFinder to find the cell and the barycentric
coordinates of x and calls evaluate_bary!.)
"""
function evaluate!(
	result,
	PE::PointEvaluator,
	x;
	kwargs...
)
	# find correct cell (start from cell of last evaluation)
	item = gFindLocal!(PE.xref, PE.CF, x; icellstart = PE.lastitem, kwargs...)
	@assert item > 0
	PE.lastitem = item

	## find cell geometry id
	j = PE.eval_selector(item)

	## evaluate
	PE.evaluator_bary(result, PE.BE_args[j], PE.L2G[j], PE.xref, item)
end



"""
````
function eval_func_bary(PE::PointEvaluator)
````

Yields the function (result, xref, item) -> evaluate_bary!(result,PE,xref,item).
"""
function eval_func_bary(PE::PointEvaluator)
	return (result, xref, item) -> evaluate_bary!(result, PE, xref, item)
end

"""
````
function eval_func(PE::PointEvaluator)
````

Yields the function (result, x) -> evaluate!(result,PE,x).
"""
function eval_func(PE::PointEvaluator; kwargs...)
	return (result, x) -> evaluate!(result, PE, x; kwargs...)
end
