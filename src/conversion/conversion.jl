# TODO: implement, test, document
# -> return one mesh, or one mesh per type of cell? (One might be interested in treating different types differently, e.g. embedded cells)
#       |-> I am not sure if Makie can treat different elements in the same mesh...
# -> treat quadratic elements differently?
# -> check how it is done in FerriteViz -> will this be needed at all?

include("typeconversions.jl")

"""
    convert_to_makie_mesh(grid; kwargs...)

Return a mesh that can be used for plotting with Makie.jl. 

Keyword arguments:
 - `displ`: a `Vector{Vec{dim, T}}` containing the nodal displacements to plot a deformed grid. The default is undeformed.
 - `cellset`: to plot a subdomain of the grid, the default is the entire `grid`.
 - `disconnectcells`: to signify if the mesh will be used to plot a discontinous field. The default is `false`.

"""
function convert_to_makie_mesh(
    grid::Grid{dim}; 
    displ::Union{Vector{Vec{dim, T}}, Nothing} = nothing, 
    cellset::Union{Vector{C},OrderedSet{Int},String}=getcells(grid),
    disconnectcells::Bool=false
    ) where {dim, T, C<:Ferrite.AbstractCell}
    celltype = typeof(first(grid.cells))
    celltype == Ferrite.Line || celltype == Ferrite.QuadraticLine ? throw("Makie.jl currently only supports 2D and 3D points") : nothing; # TODO: possible work around in the future?
    return _convert_to_makie_mesh(grid, displ, disconnectcells, cellset)
end

function _convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{Vec{dim, T}}, Nothing}, disconnectcells::Bool, cellset::OrderedSet{Int}) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, getcells(grid)[collect(cellset)])
end

function _convert_to_makie_mesh(grid::Grid{dim}, displ::Union{Vector{Vec{dim, T}}, Nothing}, disconnectcells::Bool, cellset::String) where {dim, T}
    return _convert_to_makie_mesh(grid, displ, disconnectcells, getcellset(grid, cellset))
end

function _convert_to_makie_mesh(grid::Grid{dim}, displ::Vector{Vec{dim, T}}, disconnectcells::Bool, cellset::Vector{<:Ferrite.AbstractCell}) where {dim, T}
    @assert length(displ) == getnnodes(grid)
    nodes = [(n.x + u) for (n, u) in zip(grid.nodes, displ)]
    return _convert_to_makie_mesh(nodes, Val(disconnectcells), cellset)
end

function _convert_to_makie_mesh(grid::Grid{dim}, ::Nothing, disconnectcells::Bool, cellset::Vector{<:Ferrite.AbstractCell}) where {dim}
    nodes = [n.x for n in grid.nodes]
    return _convert_to_makie_mesh(nodes, Val(disconnectcells), cellset)
end

function _convert_cells_to_makie(grid::Grid{dim}, cellset::Vector{Int}) where {dim}
    return vcat( _convert_cell_to_makie.(get_cells(grid, cellset))...) # TODO: Optimize: Preallocate final vector and write to it (in a loop?)?
end

function _convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, ::Val{false}, cells::Vector{C}) where {dim, T, C<:Ferrite.AbstractCell}
    makienodes = _convert_vec_to_makie.(nodes)
    makiefacepercell, CM = _get_makie_type_data(C)
    makiecells = Vector{CM}(undef, length(cells) * makiefacepercell)
    for (i, cell) in enumerate(cells)
        makiecells[1+(i-1)*makiefacepercell:i*makiefacepercell] .= _convert_cell_to_makie(cell)
    end
    return Makie.GeometryBasics.Mesh(makienodes, makiecells)
end

function _convert_to_makie_mesh(nodes::Vector{Vec{dim, T}}, ::Val{true}, cells::Vector{C}) where {dim, T, C<:Ferrite.AbstractCell}
    makiefacepercell, CM = _get_makie_type_data(C)
    nodespercell = Ferrite.nnodes(first(cells))
    makienodes = Vector{Makie.GeometryBasics.Point{dim,Float64}}(undef, length(cells)*nodespercell)
    makiecells = Vector{CM}(undef, length(cells) * makiefacepercell)
    for (i, cell) in enumerate(cells)
        noderange = 1+(i-1)*nodespercell:i*nodespercell
        makienodes[noderange] .= _convert_vec_to_makie.(nodes[[cell.nodes...]]) # Copy nodes (coordinates) of the current cell into ``makienodes``
        disconnected_cell = C(Tuple(n for n in noderange)) # Create a new cell using the copied nodes
        makiecells[1+(i-1)*makiefacepercell:i*makiefacepercell] .= _convert_cell_to_makie(disconnected_cell) # Convert the new cell to Makie
    end
    return Makie.GeometryBasics.Mesh(makienodes, makiecells)
end

"""
    disconnect_field(field, grid; cellset)

Returns a vector of values that can be passed to a `Makie.jl` plotting function to plot a discontinous (across element boundaries) field.

Keyword arguments:
 - `cellset`: to plot a subdomain of the grid, the default is the entire `grid`.

`mesh` should be a mesh that was generated using the `convert_to_makie_mesh(grid; disconnectcells=true)` function.
For constant values per cell `field` should be a `Vector{Float64}` with the length `Ferrite.getncells(grid)`. 
Otherwise a `Matrix` can be passed with rows corresponding the the `cellid::Int` and the columns to the value of the nodes 
corresponding to that cell.
Example matrix:
```julia
field = [
    0.44 0.29 # ... # data for nodes 1, 2 ... of element 1
    0.98 0.48 # ... # data for nodes 1, 2 ... of element 2
    0.32 0.55 # ... # data for nodes 1, 2 ... of element 3
]
```
"""
function disconnect_field(
    field::Union{Vector{Float64}, Matrix{Float64}},
    grid::Grid{dim};
    cellset::Union{OrderedSet{Int},String} = OrderedSet{Int}(1:Ferrite.getncells(grid))
    ) where {dim}
    return _disconnect_field(field, grid, cellset)
end

function _disconnect_field(field::Union{Vector{Float64}, Matrix{Float64}}, grid::Grid{dim}, cellset::String) where {dim}
    return _disconnect_field(field, grid, Ferrite.getcellset(grid, cellset))
end

function _disconnect_field(field::Union{Vector{Float64}, Matrix{Float64}}, grid::Grid{dim}, cellset::OrderedSet{Int}) where {dim}
    allcells = (Ferrite.getncells(grid) == length(cellset))
    return _disconnect_field(field, grid, cellset, Val(allcells))
end

# TODO: Benchmark for ... 
function _disconnect_field(field::Union{Vector{Float64}, Matrix{Float64}}, grid::Grid{dim}, cellset::OrderedSet{Int}, allcells::V) where {dim, V<:Union{Val{true},Val{false}}}
    _field = @view field[collect(cellset),:]
    nodespercell = Ferrite.nnodes(first(getcells(grid)))
    field_disconnected = zeros(length(cellset)*nodespercell)
    for (i, cellid) in enumerate(cellset)
        noderange = 1+(i-1)*nodespercell:i*nodespercell
        ind = _get_field_index_for_disconnection(i, cellid, allcells)
        field_disconnected[noderange] .= _field[ind,:] # Keep DOF order consistent with _convert_to_makie_mesh(nodes, ::Val{true}, cells)
    end
    return field_disconnected
end

_get_field_index_for_disconnection(i::Int, _::Int, ::Val{false}) = i
_get_field_index_for_disconnection(_::Int, cellid::Int, ::Val{true}) = cellid