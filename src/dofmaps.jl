"""
$(TYPEDEF)

Dofmaps are stored as an ExtendableGrids.AbstractGridAdjacency in the finite element space and collect
information with respect to different AssemblyTypes. They are generated automatically on demand and the dofmaps
associated to each subtype can be accessed via FESpace[DofMap].
"""
abstract type DofMap <: AbstractGridAdjacency end
abstract type CellDofs <: DofMap end
abstract type FaceDofs <: DofMap end
abstract type EdgeDofs <: DofMap end
abstract type BFaceDofs <: DofMap end
abstract type BEdgeDofs <: DofMap end
abstract type CellDofsParent <: DofMap end
abstract type FaceDofsParent <: DofMap end
abstract type EdgeDofsParent <: DofMap end
abstract type BFaceDofsParent <: DofMap end
abstract type BEdgeDofsParent <: DofMap end

const DofMapTypes{Ti} = Union{VariableTargetAdjacency{Ti}, SerialVariableTargetAdjacency{Ti}, Array{Ti, 2}}

UCG4DofMap(::Type{CellDofs}) = UniqueCellGeometries
UCG4DofMap(::Type{FaceDofs}) = UniqueFaceGeometries
UCG4DofMap(::Type{EdgeDofs}) = UniqueEdgeGeometries
UCG4DofMap(::Type{BFaceDofs}) = UniqueBFaceGeometries
UCG4DofMap(::Type{BEdgeDofs}) = UniqueBEdgeGeometries

SuperItemNodes4DofMap(::Type{CellDofs}) = CellNodes
SuperItemNodes4DofMap(::Type{FaceDofs}) = FaceNodes
SuperItemNodes4DofMap(::Type{EdgeDofs}) = EdgeNodes
SuperItemNodes4DofMap(::Type{BFaceDofs}) = FaceNodes
SuperItemNodes4DofMap(::Type{BEdgeDofs}) = EdgeNodes

ItemGeometries4DofMap(::Type{CellDofs}) = CellGeometries
ItemGeometries4DofMap(::Type{FaceDofs}) = FaceGeometries
ItemGeometries4DofMap(::Type{EdgeDofs}) = EdgeGeometries
ItemGeometries4DofMap(::Type{BFaceDofs}) = BFaceGeometries
ItemGeometries4DofMap(::Type{BEdgeDofs}) = BEdgeGeometries

ItemEdges4DofMap(::Type{CellDofs}) = CellEdges
ItemEdges4DofMap(::Type{FaceDofs}) = FaceEdges
ItemEdges4DofMap(::Type{BFaceDofs}) = FaceEdges

Sub2Sup4DofMap(::Type{<:DofMap}) = nothing
Sub2Sup4DofMap(::Type{BFaceDofs}) = BFaceFaces
Sub2Sup4DofMap(::Type{BEdgeDofs}) = BEdgeEdges

Dofmap4AssemblyType(::Type{ON_CELLS}) = CellDofs
Dofmap4AssemblyType(::Type{<:ON_FACES}) = FaceDofs
Dofmap4AssemblyType(::Type{ON_BFACES}) = BFaceDofs
Dofmap4AssemblyType(::Type{<:ON_EDGES}) = EdgeDofs
Dofmap4AssemblyType(::Type{ON_BEDGES}) = BEdgeDofs

ParentDofmap4Dofmap(::Type{CellDofs}) = CellDofsParent
ParentDofmap4Dofmap(::Type{FaceDofs}) = FaceDofsParent
ParentDofmap4Dofmap(::Type{EdgeDofs}) = EdgeDofsParent
ParentDofmap4Dofmap(::Type{BFaceDofs}) = BFaceDofsParent
ParentDofmap4Dofmap(::Type{BEdgeDofs}) = BEdgeDofsParent

EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{ON_CELLS}) = ON_CELLS
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{<:ON_FACES}) = ON_FACES
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{ON_BFACES}) = ON_BFACES
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{<:ON_EDGES}) = ON_EDGES
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{ON_BEDGES}) = ON_BEDGES

EffAT4AssemblyType(::Type{ON_FACES}, ::Type{ON_CELLS}) = nothing
EffAT4AssemblyType(::Type{ON_FACES}, ::Type{<:ON_FACES}) = ON_CELLS
EffAT4AssemblyType(::Type{ON_FACES}, ::Type{<:ON_EDGES}) = ON_FACES
EffAT4AssemblyType(::Type{ON_FACES}, ::Type{<:ON_BEDGES}) = ON_BFACES

EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{ON_CELLS}) = nothing
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{ON_BFACES}) = ON_CELLS
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{<:ON_FACES}) = nothing
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{<:ON_EDGES}) = nothing
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{<:ON_BEDGES}) = nothing

EffAT4AssemblyType(::Type{ON_EDGES}, ::Type{ON_CELLS}) = nothing
EffAT4AssemblyType(::Type{ON_EDGES}, ::Type{<:ON_FACES}) = nothing
EffAT4AssemblyType(::Type{ON_EDGES}, ::Type{<:ON_EDGES}) = ON_CELLS


abstract type DofType end
abstract type DofTypeNode <: DofType end        # parsed from 'N' or 'n' (nodal continuous dof)
abstract type DofTypeFace <: DofType end        # parsed from 'F' or 'f' (face continuous dof)
abstract type DofTypeEdge <: DofType end        # parsed from 'E' or 'e' (edge continuous dof)
abstract type DofTypeInterior <: DofType end    # parsed from 'I' or 'i' (interior dof)
abstract type DofTypePCell <: DofType end       # parsed from 'C' or 'c' (parent cell dof, only needed by P0 element for BFaceDofs)

function DofType(c::Char)
	if lowercase(c) == 'n'
		return DofTypeNode
	elseif lowercase(c) == 'f'
		return DofTypeFace
	elseif lowercase(c) == 'e'
		return DofTypeEdge
	elseif lowercase(c) == 'i'
		return DofTypeInterior
	elseif lowercase(c) == 'c'
		return DofTypePCell
	else
		@error "No DofType available to parse from $(lowercase(c))"
	end
end

const dofmap_type_chars = ['E', 'N', 'I', 'C', 'F', 'e', 'n', 'i', 'c', 'f']
const dofmap_number_chars = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']

struct DofMapPatternSegment
	type::Type{<:DofType}
	each_component::Bool    # single dof for all components or one for each component ?
	ndofs::Int              # how many dofs 
end

# splits pattern into pairs of single chars and Ints
function parse_pattern(pattern::String)
	pairs = Array{DofMapPatternSegment, 1}([])
	for j ∈ 1:length(pattern)
		if pattern[j] in dofmap_type_chars
			k = j + 1
			while (k < length(pattern))
				if (pattern[k+1] in dofmap_number_chars)
					k += 1
				else
					break
				end
			end
			# @show pattern[j] SubString(pattern,j+1,k)
			push!(pairs, DofMapPatternSegment(DofType(pattern[j]), isuppercase(pattern[j]), tryparse(Int, SubString(pattern, j + 1, k))))
		end
	end
	return pairs
end


struct ParsedDofMap
	segments::Array{DofMapPatternSegment, 1}
	ndofs_node4c::Int
	ndofs_face4c::Int
	ndofs_edge4c::Int
	ndofs_interior4c::Int
	ndofs_pcell4c::Int
	ndofs_node::Int
	ndofs_face::Int
	ndofs_edge::Int
	ndofs_interior::Int
	ndofs_pcell::Int
	ndofs_total::Int
end

get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeNode}) = P.ndofs_node4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeEdge}) = P.ndofs_edge4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeFace}) = P.ndofs_face4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeInterior}) = P.ndofs_interior4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypePCell}) = P.ndofs_pcell4c
get_ndofs(P::ParsedDofMap, ::Type{DofTypeNode}) = P.ndofs_node
get_ndofs(P::ParsedDofMap, ::Type{DofTypeEdge}) = P.ndofs_edge
get_ndofs(P::ParsedDofMap, ::Type{DofTypeFace}) = P.ndofs_face
get_ndofs(P::ParsedDofMap, ::Type{DofTypeInterior}) = P.ndofs_interior
get_ndofs(P::ParsedDofMap, ::Type{DofTypePCell}) = P.ndofs_pcell
get_ndofs(P::ParsedDofMap) = P.ndofs_total

function ParsedDofMap(pattern::String, ncomponents, EG::Type{<:AbstractElementGeometry})
	segments = parse_pattern(pattern)
	ndofs_node::Int = 0
	ndofs_face::Int = 0
	ndofs_edge::Int = 0
	ndofs_interior::Int = 0
	ndofs_pcell::Int = 0
	ndofs_node4c::Int = 0
	ndofs_face4c::Int = 0
	ndofs_edge4c::Int = 0
	ndofs_interior4c::Int = 0
	ndofs_pcell4c::Int = 0
	for j ∈ 1:length(segments)
		if segments[j].type <: DofTypeNode
			ndofs_node += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
			if segments[j].each_component
				ndofs_node4c += segments[j].ndofs
			end
		elseif segments[j].type <: DofTypeFace
			ndofs_face += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
			if segments[j].each_component
				ndofs_face4c += segments[j].ndofs
			end
		elseif segments[j].type <: DofTypeEdge
			ndofs_edge += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
			if segments[j].each_component
				ndofs_edge4c += segments[j].ndofs
			end
		elseif segments[j].type <: DofTypeInterior
			ndofs_interior += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
			if segments[j].each_component
				ndofs_interior4c += segments[j].ndofs
			end
		elseif segments[j].type <: DofTypePCell
			ndofs_pcell += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
			if segments[j].each_component
				ndofs_pcell4c += segments[j].ndofs
			end
		end
	end
	ndofs_total = num_nodes(EG) * ndofs_node + num_faces(EG) * ndofs_face + num_edges(EG) * ndofs_edge + ndofs_interior
	return ParsedDofMap(segments, ndofs_node4c, ndofs_face4c, ndofs_edge4c, ndofs_interior4c, ndofs_pcell4c, ndofs_node, ndofs_face, ndofs_edge, ndofs_interior, ndofs_pcell, ndofs_total)
end


function Dofmap4AssemblyType(FES::FESpace, AT::Type{<:AssemblyType})
	return FES[Dofmap4AssemblyType(EffAT4AssemblyType(assemblytype(FES), AT))]
end

function init_dofmap_from_pattern!(FES::FESpace{Tv, Ti, FEType, APT}, DM::Type{<:DofMap}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
	## Beware: Automatic broken DofMap generation currently only reliable for CellDofs

	@debug "Generating $DM for $(FES.name)"
	## prepare dofmap patterns
	xgrid = FES.dofgrid
	EG = xgrid[UCG4DofMap(DM)]
	ncomponents::Int = get_ncomponents(FEType)
	need_faces::Bool = false
	need_edges::Bool = false
	maxdofs4item::Int = 0
	dofmap4EG::Array{ParsedDofMap, 1} = Array{ParsedDofMap, 1}(undef, length(EG))
	has_interiordofs = zeros(Bool, length(EG))
	for j ∈ 1:length(EG)
		pattern = get_dofmap_pattern(FEType, DM, EG[j])
		dofmap4EG[j] = ParsedDofMap(pattern, ncomponents, EG[j])
		maxdofs4item = max(maxdofs4item, get_ndofs(dofmap4EG[j]))
		if get_ndofs(dofmap4EG[j], DofTypeFace) > 0
			need_faces = true
		end
		if get_ndofs(dofmap4EG[j], DofTypeEdge) > 0
			need_edges = true
		end
		if get_ndofs(dofmap4EG[j], DofTypeInterior) > 0
			has_interiordofs[j] = true
		end
	end

	## prepare data
	nnodes::Int = size(xgrid[Coordinates], 2)
	xItemNodes::Adjacency{Ti} = xgrid[SuperItemNodes4DofMap(DM)]
	xItemGeometries::GridEGTypes = xgrid[ItemGeometries4DofMap(DM)]
	if need_faces
		xItemFaces::Adjacency{Ti} = xgrid[CellFaces]
		nfaces = num_sources(xgrid[FaceNodes])
	end
	if need_edges
		xItemEdges::Adjacency{Ti} = xgrid[ItemEdges4DofMap(DM)]
		nedges = num_sources(xgrid[EdgeNodes])
	end
	nitems::Int = num_sources(xItemNodes)
	offset4component = FES.coffset

	## generate dofmap from patterns
	sub2sup = nothing
	nsubitems::Int = 0
	if Sub2Sup4DofMap(DM) !== nothing
		Sub2SupItems = xgrid[Sub2Sup4DofMap(DM)]
		sub2sup = (x) -> Sub2SupItems[x]
		nsubitems = length(Sub2SupItems)
	else
		sub2sup = (x) -> x
		nsubitems = nitems
	end

	nnodes4EG::Array{Int, 1} = zeros(Int, length(EG))
	nfaces4EG::Array{Int, 1} = zeros(Int, length(EG))
	nedges4EG::Array{Int, 1} = zeros(Int, length(EG))
	for j ∈ 1:length(EG)
		nnodes4EG[j] = num_nodes(EG[j])
		nfaces4EG[j] = num_faces(EG[j])
		nedges4EG[j] = num_edges(EG[j])
	end

	itemEG = EG[1]
	iEG::Int = 1

	## get total dofnumber for pre--allocation
	dofmap_totallength = 0
	colstarts = zeros(Ti, nsubitems + 1)
	for subitem ∈ 1:nsubitems
		itemEG = xItemGeometries[subitem]
		if length(EG) > 1
			iEG = findfirst(isequal(itemEG), EG)
		end
		colstarts[subitem] = dofmap_totallength + 1
		dofmap_totallength += get_ndofs(dofmap4EG[iEG])
	end
	colstarts[end] = dofmap_totallength + 1

	if FES.broken
		FES[DM] = SerialVariableTargetAdjacency(colstarts)
		return FES[DM]
	end

	colentries = zeros(Ti, dofmap_totallength)
	xItemDofs = VariableTargetAdjacency{Ti}(colentries, colstarts)

	cpattern::Array{DofMapPatternSegment, 1} = dofmap4EG[1].segments
	l::Int = 0
	offset::Int = 0
	pos::Int = 0
	q::Int = 0
	item_with_interiordofs::Int = 0
	for subitem ∈ 1:nsubitems
		itemEG = xItemGeometries[subitem]
		item = sub2sup(subitem)
		if length(EG) > 1
			iEG = findfirst(isequal(itemEG), EG)
		end
		cpattern = dofmap4EG[iEG].segments
		l = length(cpattern)
		if has_interiordofs[iEG]
			item_with_interiordofs += 1
		end
		for c ∈ 1:ncomponents
			offset = (c - 1) * offset4component
			for k ∈ 1:l
				q = cpattern[k].ndofs
				if cpattern[k].type <: DofTypeNode && cpattern[k].each_component
					for n ∈ 1:nnodes4EG[iEG]
						for m ∈ 1:q
							pos += 1
							xItemDofs.colentries[pos] = xItemNodes[n, item] + offset + (m - 1) * nnodes
						end
					end
					offset += nnodes * q
				elseif cpattern[k].type <: DofTypeFace && cpattern[k].each_component
					for n ∈ 1:nfaces4EG[iEG]
						for m ∈ 1:q
							pos += 1
							xItemDofs.colentries[pos] = xItemFaces[n, item] + offset + (m - 1) * nfaces
						end
					end
					offset += nfaces * q
				elseif cpattern[k].type <: DofTypeEdge && cpattern[k].each_component
					for n ∈ 1:nedges4EG[iEG]
						for m ∈ 1:q
							pos += 1
							xItemDofs.colentries[pos] = xItemEdges[n, item] + offset + (m - 1) * nedges
						end
					end
					offset += nedges * q
				elseif cpattern[k].type <: DofTypeInterior && cpattern[k].each_component
					for m ∈ 1:q
						pos += 1
						xItemDofs.colentries[pos] = sub2sup(item_with_interiordofs) + offset
						offset += nitems
					end
				end
			end
		end
		offset = ncomponents * offset4component
		for k ∈ 1:l
			q = cpattern[k].ndofs
			if cpattern[k].type <: DofTypeFace && !cpattern[k].each_component
				for n ∈ 1:nfaces4EG[iEG]
					for m ∈ 1:q
						pos += 1
						xItemDofs.colentries[pos] = xItemFaces[n, item] + offset + (m - 1) * nfaces
					end
				end
				offset += nfaces * q
			elseif cpattern[k].type <: DofTypeEdge && !cpattern[k].each_component
				for n ∈ 1:nedges4EG[iEG]
					for m ∈ 1:q
						pos += 1
						xItemDofs.colentries[pos] = xItemEdges[n, item] + offset + (m - 1) * nedges
					end
				end
				offset += nedges * q
			elseif cpattern[k].type <: DofTypeInterior && !cpattern[k].each_component
				for m ∈ 1:q
					pos += 1
					xItemDofs.colentries[pos] = sub2sup(item_with_interiordofs) + offset
					offset += nitems
				end
			end
		end
	end

	FES[DM] = xItemDofs

	if FES.dofgrid !== FES.xgrid
		## assume parent relation between xgrid and dofgrid
		@assert FES.dofgrid[ParentGrid] == FES.xgrid "xgrid is not the parent grid of dofgrid !!!"
		@assert FES.dofgrid[ParentGridRelation] == SubGrid "dofgrid is not a subgrid of xgrid !!!"
		## lift subgrid dofmap to parent xgrid to allow assembly on common parent grid
		## by constructing a variable target adjacency with empty dofs for parent items in xgrid
		## that are not in the dofgrid --> this allows to use the dofmap in assembly over full xgrid
		xItemDofsLifted = VariableTargetAdjacency(Ti)
		if DM == CellDofs
			parentitems = FES.dofgrid[CellParents]
		elseif DM == FaceDofs
			parentitems = FES.dofgrid[FaceParents]
		elseif DM == EdgeDofs
			parentitems = FES.dofgrid[EdgeParents]
		elseif DM == BFaceDofs
			parentitems = FES.dofgrid[BFaceParents]
		elseif DM == BEdgeDofs
			parentitems = FES.dofgrid[BEdgeParents]
		end
		DMP = ParentDofmap4Dofmap(DM)
		subitem = 1
		for i = 1 : num_sources(xItemDofs)
			t = parentitems[i]
			for j = subitem : t-1
				append!(xItemDofsLifted, [])
				subitem += 1
			end
			append!(xItemDofsLifted, view(xItemDofs, :, i))
			subitem += 1
		end

		FES[DMP] = xItemDofsLifted
	end
	return FES[DM]
end


function init_broken_dofmap!(FES::FESpace{Tv, Ti, FEType, APT}, DM::Union{Type{BFaceDofs}, Type{FaceDofs}}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}

	## prepare dofmap patterns
	xgrid = FES.dofgrid
	cellEG = xgrid[UniqueCellGeometries]
	EG = xgrid[UCG4DofMap(DM)]
	ncomponents::Int = get_ncomponents(FEType)
	need_nodes = false
	need_faces = false
	need_edges = false
	maxdofs4item::Int = 0
	dofmap4cellEG::Array{ParsedDofMap, 1} = Array{ParsedDofMap, 1}(undef, length(cellEG))
	dofmap4EG::Array{ParsedDofMap, 1} = Array{ParsedDofMap, 1}(undef, length(EG))
	for j ∈ 1:length(cellEG)
		pattern = get_dofmap_pattern(FEType, CellDofs, cellEG[j])
		dofmap4cellEG[j] = ParsedDofMap(pattern, ncomponents, cellEG[j])
	end
	for j ∈ 1:length(EG)
		pattern = get_dofmap_pattern(FEType, DM, EG[j])
		dofmap4EG[j] = ParsedDofMap(pattern, ncomponents, EG[j])
		maxdofs4item = max(maxdofs4item, get_ndofs(dofmap4EG[j]))
		if get_ndofs(dofmap4EG[j], DofTypeNode) > 0
			need_nodes = true
		end
		if get_ndofs(dofmap4EG[j], DofTypeInterior) > 0
			need_faces = true
		end
		if get_ndofs(dofmap4EG[j], DofTypeEdge) > 0
			need_edges = true
		end
	end

	xFaceCells = xgrid[FaceCells]
	xCellNodes = xgrid[CellNodes]
	xCellDofs = FES[CellDofs]
	xFaceNodes = xgrid[SuperItemNodes4DofMap(DM)]
	xCellGeometries = xgrid[CellGeometries]
	xFaceGeometries = xgrid[ItemGeometries4DofMap(DM)]
	xFaceDofs = VariableTargetAdjacency(Ti)
	if DM == BFaceDofs
		xRealFace = xgrid[BFaceFaces]
		nfaces = length(xRealFace)
	elseif DM == FaceDofs
		nfaces = num_sources(xFaceNodes)
		xRealFace = 1:nfaces
	end

	if need_faces
		xCellFaces = xgrid[CellFaces]
	end
	if need_edges
		xFaceEdges = xgrid[FaceEdges]
		xCellEdges = xgrid[CellEdges]
		local_edges = zeros(Int, max_num_targets_per_source(xFaceEdges))
	end

	cpattern::Array{DofMapPatternSegment, 1} = dofmap4EG[1].segments
	ccellpattern::Array{DofMapPatternSegment, 1} = dofmap4cellEG[1].segments
	iEG::Int = 1
	local_nodes = zeros(Int, max_num_targets_per_source(xFaceNodes))
	local_face::Int = 0
	local_dofs = zeros(Int, 2 * max_num_targets_per_source(xCellDofs))
	nldofs::Int = 0
	celldofoffset::Int = 0
	node::Int = 0
	face::Int = 0
	cell::Int = 0
	nfacenodes::Int = 0
	nfaceedges::Int = 0
	ncellnodes::Int = 0
	ncelledges::Int = 0
	ncellfaces::Int = 0
	pos::Int = 0
	cEG = xCellGeometries[1]
	for f ∈ 1:nfaces
		face = xRealFace[f]
		cEG = xFaceGeometries[f]
		iEG = findfirst(isequal(cEG), EG)
		cpattern = dofmap4EG[iEG].segments
		nldofs = 0
		for icell ∈ 1:2
			cell = xFaceCells[icell, face]
			if cell > 0
				cEG = xCellGeometries[cell]
				iEG = findfirst(isequal(cEG), cellEG)
				ccellpattern = dofmap4cellEG[iEG].segments

				## get local nodes/faces/edges dofs for global face f
				if need_nodes
					ncellnodes = num_targets(xCellNodes, cell)
					nfacenodes = num_targets(xFaceNodes, face)
					for k ∈ 1:nfacenodes
						node = xFaceNodes[k, face]
						pos = 1
						while xCellNodes[pos, cell] != node
							pos += 1
						end
						local_nodes[k] = pos
					end
				end
				if need_faces
					ncellfaces = num_targets(xCellFaces, cell)
					pos = 1
					while xCellFaces[pos, cell] != face
						pos += 1
					end
					local_face = pos
				end
				if need_edges
					ncelledges = num_targets(xCellEdges, cell)
					nfaceedges = num_targets(xFaceEdges, face)
					for k ∈ 1:nfaceedges
						edge = xFaceEdges[k, face]
						pos = 1
						while xCellEdges[pos, cell] != edge
							pos += 1
						end
						local_edges[k] = pos
					end
				end

				## get face dofs of cell and write them into local_dofs
				## assuming that CellDofs and FaceDofs patterns are consistent (e.g. "N1F1" and "N1I1")

				## get each-component dofs on cell
				celldofoffset = 0
				for c ∈ 1:ncomponents
					for s ∈ 1:length(cpattern)
						q = cpattern[s].ndofs
						if cpattern[s].type <: DofTypePCell && cpattern[s].each_component
							for dof ∈ 1:q
								nldofs += 1
								local_dofs[nldofs] = xCellDofs[celldofoffset+dof, cell]
							end
							celldofoffset += q
						elseif cpattern[s].type <: DofTypeNode && cpattern[s].each_component
							for k ∈ 1:nfacenodes
								for dof ∈ 1:q
									nldofs += 1
									local_dofs[nldofs] = xCellDofs[celldofoffset+(local_nodes[k]-1)*q+dof, cell]
								end
							end
							celldofoffset += q * ncellnodes
						elseif cpattern[s].type <: DofTypeInterior && cpattern[s].each_component
							for dof ∈ 1:q
								nldofs += 1
								local_dofs[nldofs] = xCellDofs[celldofoffset+(local_face-1)*q+dof, cell]
							end
							celldofoffset += q * ncellfaces
						elseif cpattern[s].type <: DofTypeEdge && cpattern[s].each_component
							for k ∈ 1:nfaceedges
								for dof ∈ 1:q
									nldofs += 1
									local_dofs[nldofs] = xCellDofs[celldofoffset+(local_edges[k]-1)*q+dof, cell]
								end
							end
							celldofoffset += q * ncelledges
						end
					end
					# increase celldofoffset for interior component dofs in cell pattern
					for s ∈ 1:length(ccellpattern)
						if ccellpattern[s].type <: DofTypeInterior && ccellpattern[s].each_component
							celldofoffset += ccellpattern[s].ndofs
						end
					end
				end

				## add single component dofs on cell
				for s ∈ 1:length(cpattern)
					q = cpattern[s].ndofs
					if cpattern[s].type <: DofTypeInterior && !cpattern[s].each_component
						for dof ∈ 1:q
							nldofs += 1
							local_dofs[nldofs] = xCellDofs[celldofoffset+(local_face-1)*q+dof, cell]
						end
						celldofoffset += q
					elseif cpattern[s].type <: DofTypePCell && !cpattern[s].each_component
						for dof ∈ 1:q
							nldofs += 1
							local_dofs[nldofs] = xCellDofs[celldofoffset+dof, cell]
						end
						celldofoffset += q
					elseif cpattern[s].type <: DofTypeNode && !cpattern[s].each_component
						for k ∈ 1:nfacenodes
							for dof ∈ 1:q
								nldofs += 1
								local_dofs[nldofs] = xCellDofs[celldofoffset+(local_nodes[k]-1)*q+dof, cell]
							end
						end
						celldofoffset += q * ncellnodes
					elseif cpattern[s].type <: DofTypeInterior && !cpattern[s].each_component
						for dof ∈ 1:q
							nldofs += 1
							local_dofs[nldofs] = xCellDofs[celldofoffset+(local_face-1)*q+dof, cell]
						end
						celldofoffset += q
					elseif cpattern[s].type <: DofTypeEdge && !cpattern[s].each_component
						for k ∈ 1:nfaceedges
							for dof ∈ 1:q
								nldofs += 1
								local_dofs[nldofs] = xCellDofs[celldofoffset+(local_edges[k]-1)*q+dof, cell]
							end
						end
						celldofoffset += q * ncelledges
					end
				end
			end
		end
		append!(xFaceDofs, local_dofs[1:nldofs])
	end
	FES[DM] = xFaceDofs
end


function init_dofmap!(FES::FESpace, DM::Type{<:DofMap})
	@debug "Generating dofmap $DM for FESpace $(FES.name)"

	if (FES.broken == true) && (DM != CellDofs)
		## dofmap needs to include all (e.g. face) dofs from neighbouring cells
		init_broken_dofmap!(FES, DM)
	else
		init_dofmap_from_pattern!(FES, DM)
	end
end

