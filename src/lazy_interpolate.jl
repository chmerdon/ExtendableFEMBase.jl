

"""
````
function lazy_interpolate!(
	target::FEVectorBlock{T1,Tv,Ti},
	source::FEVectorBlock{T2,Tv,Ti};
	operator = Identity,
	postprocess = nothing,
	xtrafo = nothing,
	items = [],
	not_in_domain_value = 1e30,
	use_cellparents::Bool = false) where {T1,T2,Tv,Ti}
````

Interpolates (operator-evaluations of) the given finite element function into the finite element space assigned to the target FEVectorBlock. 
(Currently not the most efficient way as it is based on the PointEvaluation pattern and cell search. If CellParents
are available in the grid components of the target grid, these parent cell information can be used to improve the
search. To activate this put 'use_cellparents' = true). By some given kernel function (result, input, qpinfo)
the operator evaluation (=input) can be further postprocessed.

Note: discontinuous quantities at vertices of the target grid will be evaluted in the first found cell of the
source grid. No averaging is performed.
"""
function lazy_interpolate!(
	target::FEVectorBlock{T1, Tv, Ti},
	source,
	operators = [(1, Identity)];
	postprocess = standard_kernel,
	xtrafo = nothing,
	items = [],
	resultdim = get_ncomponents(eltype(target.FES)),
	not_in_domain_value = 1e30,
	start_cell = 1,
	only_localsearch = false,
	use_cellparents::Bool = false,
	eps = 1e-13,
	kwargs...) where {T1, Tv, Ti}

	# wrap point evaluation into function that is put into normal interpolate!
	xgrid = source[1].FES.xgrid
	xdim_source::Int = size(xgrid[Coordinates], 1)
	xdim_target::Int = size(target.FES.xgrid[Coordinates], 1)
	if xdim_source != xdim_target
		@assert xtrafo !== nothing "grids have different coordinate dimensions, need xtrafo!"
	end
	PE = PointEvaluator(postprocess, operators, source)
	xref = zeros(Tv, xdim_source)
	x_source = zeros(Tv, xdim_source)
	cell::Int = start_cell
	lastnonzerocell::Int = start_cell
	same_cells::Bool = xgrid == target.FES.xgrid
	CF::CellFinder{Tv, Ti} = CellFinder(xgrid)

	if same_cells || use_cellparents == true
		if same_cells
			xCellParents = 1:num_cells(target.FES.xgrid)
		else
			xCellParents::Array{Ti, 1} = target.FES.xgrid[CellParents]
		end
		function point_evaluation_parentgrid!(result, qpinfo)
			x = qpinfo.x
			cell = xCellParents[qpinfo.cell]
			if xtrafo !== nothing
				xtrafo(x_source, x)
				cell = gFindLocal!(xref, CF, x_source; icellstart = cell, eps = eps, trybrute = !only_localsearch)
			else
				cell = gFindLocal!(xref, CF, x; icellstart = cell, eps = eps, trybrute = !only_localsearch)
			end
			evaluate!(result, PE, xref, cell)
			return nothing
		end
		fe_function = point_evaluation_parentgrid!
	else
		function point_evaluation_arbitrarygrids!(result, qpinfo)
			x = qpinfo.x
			if xtrafo !== nothing
				xtrafo(x_source, x)
				cell = gFindLocal!(xref, CF, x_source; icellstart = lastnonzerocell, eps = eps, trybrute = !only_localsearch)
			else
				cell = gFindLocal!(xref, CF, x; icellstart = lastnonzerocell, eps = eps, trybrute = !only_localsearch)
			end
			if cell == 0
				fill!(result, not_in_domain_value)
			else
				evaluate_bary!(result, PE, xref, cell)
				lastnonzerocell = cell
			end
			return nothing
		end
		fe_function = point_evaluation_arbitrarygrids!
	end
	interpolate!(target, ON_CELLS, fe_function; resultdim = resultdim, items = items, kwargs...)
end
