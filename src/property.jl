"""
    get_moment(grid, order::Int; kwargs...)

Return the moment of a `grid` or a `cellset`. The `order ∈ {0, 1, 2}` correspondes to the 
volume, centre of gravity and second moment of inertia respectively. The moment is calculated using the default `QuadratureRule`
for the cells of the grid.

Keyword arguments:
 - `cellset`: to determine the moment of a subdomain of the grid, the default is the entire `grid`.
 - `refpoint`: to change the reference point, the default is the the origin. Only applicable for `order ≥ 1`.

"""
function get_moment(grid::Grid{dim}, order::Int; 
        cellset::Union{String,AbstractSet{Int}}=OrderedSet{Int}(1:getncells(grid)), 
        refpoint::Vec{dim}=zero(Vec{dim})) where {dim}
    if cellset isa String
        cellset = getcellset(grid, cellset)
    end
    cellidsbytype = _get_cellids_by_type(grid, cellset)
    return _compute_moment(grid, Val(order), cellidsbytype, refpoint)
end

function _compute_moment(grid::Grid{dim}, order::Val{ord}, cellidsbytype::Dict{DataType,OrderedSet{Int}}, refpoint::Vec{dim}) where {dim, ord}
    m = _init_moment(order, dim)
    for (T, cellids) in cellidsbytype
        cv = _init_cv(T, order)
        for cellid in cellids
            coords = getcoordinates(grid, cellid)
            reinit!(cv, coords)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, qp)
                x  = spatial_coordinate(cv, qp, coords)
                m += _compute_point_moment(x, refpoint, order) * dΩ
            end
        end
    end
    return m
end

_init_moment(::Val{0}, dim::Int) = 0.0
_init_moment(::Val{1}, dim::Int) = zero(Vec{dim})
_init_moment(::Val{2}, dim::Int) = zero(Tensor{2,dim})

function _init_cv(T::Type{<:Ferrite.AbstractCell{RefShape}}, ::Val{ord}) where {RefShape, ord}
    qr  = QuadratureRule{RefShape}(1+ord)
    ipf = Lagrange{RefShape,1}()
    ipg = geometric_interpolation(T)
    return CellValues(qr, ipf, ipg)
end

@inline _compute_point_moment(x::Vec{dim}, x̂::Vec{dim}, ::Val{0}) where {dim} = 1.0
@inline _compute_point_moment(x::Vec{dim}, x̂::Vec{dim}, ::Val{1}) where {dim} = x - x̂
@inline _compute_point_moment(x::Vec{dim}, x̂::Vec{dim}, ::Val{2}) where {dim} = (x - x̂) ⊗ (x - x̂)

###################################################################################################
###################################################################################################
# TODO: test, document
# -> Is it more efficient to not collect the coordinates in a matrix, but iterate over grid.nodes and check against current min and max?
"""
    get_coordinate_limits(grid::Grid)

Return a tuple of the minimum and maximum values for each `grid` direction.
"""
function get_coordinate_limits(grid::Grid{dim}) where {dim}
    coords = [ n.x[i] for n in grid.nodes, i in 1:dim ]
    return Tuple( (minimum(coords[:,i]), maximum(coords[:,i])) for i in 1:dim )
end

###################################################################################################
###################################################################################################
# TODO: test, document
# -> working for all dimensions like this?
# -> return two facet sets, one from the perspective of each cell set?
"""
    get_interface_between_sets(grid::Grid, set¹::AbstractSet{Int}, set²::AbstractSet{Int})

Return an `OrderedSet{FacetIndex}` of the facets building an interface between the two sets `set¹` and `set²`.
"""
function get_interface_between_sets(grid::Grid{dim}, set¹::AbstractSet{Int}, set²::AbstractSet{Int}) where {dim}
    top = ExclusiveTopology(grid)
    Γⁱⁿᵗ = OrderedSet{FacetIndex}()
    for cellid¹ in set¹
        cell¹ = grid.cells[cellid¹]
        for facet in 1:nfacets(cell¹)
            facetindex = FacetIndex(cellid¹, facet)
            neighbor = getneighborhood(top, grid, facetindex)
            isempty(neighbor) ? continue : nothing
            # length(neighbor) == 0 ? continue : nothing
            if neighbor[1][1] in set²
                push!(Γⁱⁿᵗ, facetindex)
            end
        end
    end
    return Γⁱⁿᵗ
end
function get_interface_between_sets(grid, set¹::String, set²::String)
    return get_interface_between_sets(grid, getcellset(grid, set¹), getcellset(grid, set²))
end
function get_interface_between_sets(grid, set¹::String, set²::AbstractSet{Int})
    return get_interface_between_sets(grid, getcellset(grid, set¹), set²)
end
function get_interface_between_sets(grid, set¹::AbstractSet{Int}, set²::String)
    return get_interface_between_sets(grid, set¹, getcellset(grid, set²))
end

"""
    get_dofs_from_coordinate(dh::DofHandler, x::Vec, fieldname::Symbol; radius=1e-12)

A brute force approach to finding the degrees of freedom corresponding to `fieldname` for a node at coordinate `x`. 
The node must be within a neighbourhood of radius `radius`. The default `radius` is 1e-12.
"""
function get_dofs_from_coordinate(dh::DofHandler{dim}, x::Vec{dim}, fieldname::Symbol; radius::Float64=1e-12) where {dim}
    ip = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, fieldname))
    isa(ip, ScalarInterpolation) ? dofs_per_field = 1 : dofs_per_field = dim
    nodeid, cellid, found_node = get_node_from_coordinate(dh, x; radius=radius)
    if found_node
        node_position = findfirst(x -> x == nodeid, getcells(dh.grid)[cellid].nodes)
        cell_dofs = celldofs(dh, cellid)
        node_dofs = cell_dofs[dof_range(dh, fieldname)][dofs_per_field * (node_position - 1) + 1 : dofs_per_field * node_position]
        return node_dofs
    else
        throw("No node was found with coordinate $x, try increasing the radius")
    end
end

function get_dofs_from_coordinate(dh::DofHandler{dim}, x::Vector, fieldname::Symbol; radius::Float64=1e-12) where {dim}
    return get_dofs_from_coordinate(dh, Tensors.Tensor{1,dim}(x), fieldname; radius=radius)
end

"""
    get_node_from_coordinate(dh::DofHandler, x::Vector; radius=1e-12)

Brute force search `grid` for the first node within a neighbourhood of radius `radius` around `x`.
The node id, a cell id to which the node belongs and a bool to signify a succesful search are returned.
"""
function get_node_from_coordinate(dh::DofHandler{dim}, x::Vec{dim}; radius::Float64=1e-12) where {dim}
    cells = getcells(dh.grid)
    node_position = nothing
    @inline _get_node_coords(n) = Ferrite.get_node_coordinate(dh.grid, n)
    for cellid in eachindex(cells) # possiblity for errors if cells is not indexed linearly
        node_position = findfirst(_x -> norm(_x - x) < radius, _get_node_coords.(cells[cellid].nodes))
        if !isnothing(node_position)
            nodeid = cells[cellid].nodes[node_position]
            return nodeid, cellid, true
        end
    end
    return nothing, nothing, false
end

function get_node_from_coordinate(dh::DofHandler{dim}, x::Vector; radius::Float64=1e-12) where {dim}
    return get_node_from_coordinate(dh, Tensors.Tensor{1, dim}(x); radius=radius)
end

"""
    get_coordinate_from_nodeid(id::Int, grid::Grid)

A brute force approach to find the coordinates for a node id, returns `nothing` if no id is found.
"""
function get_coordinate_from_nodeid(id::Int, grid::Grid{dim}) where {dim}
    @assert 1 <= id <= Ferrite.getnnodes(grid)
    nodeid = nothing
    for cell in Ferrite.CellIterator(grid)
        nodeid = findfirst(x -> x == id, getnodes(cell))
        if isnothing(nodeid)
            continue
        else
            return getcoordinates(cell)[nodeid]
        end
    end
    return nodeid
end

"""
!!! Limitations
    `FerriteGridUtil.get_coordinate_from_nodeid` and `FerriteGridUtil.get_dofs_from_coordinate` perform brute force
    searches, in the case that performance is key alternatives should be explored.
"""
