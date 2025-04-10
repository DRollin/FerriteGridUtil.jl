# TODO: implement, test, document
# -> to an HDF5 file?
#       |-> relatively simple, but some work to generalize for all kinds of elements etc
#       |-> add to FerriteMeshParser.jl?
# -> one of the "NASTRAN® Bulk Data" file types: .nas, .bdf, .nastran, .dat to allow import to COMSOL: https://www.comsol.com/support/learning-center/article/76161
#       |-> Not sure how these files are structured... 
#       |-> Gmsh should be able to export to a .bdf file !!! -> A function to add to Gmsh could even be used to export to other types *_*
#       |-> make this an extension of FerriteGmsh? Or just add it directly to FerriteGmsh? (so far they only translate from Gmsh to Ferrite)
# -> .stl file for whole grid and only a cellset?

"""
    save(grid::Grid{dim, C}, filepath::String) where {dim, C}

Saves a grid `grid` in an HDF5 file at `filepath`.
"""
function save(grid::Grid{dim, C}, filepath::String) where {dim, C}
    f = h5open(filepath, "w")
    _save_nodes!(f, grid)
    _save_sets!(f, grid.cellsets,  "cell sets")
    _save_sets!(f, grid.facetsets, "facet sets")
    _save_sets!(f, grid.nodesets,  "node sets")
    _save_cells!(f, grid)
    close(f)
    return grid
end

function _save_nodes!(f, grid::Grid{dim, C}) where {dim, C}
    g = create_group(f, "nodes")
    g["coords"] = [ n.x.data for n in grid.nodes ]
end

function _save_sets!(f, sets::Dict, groupname::String)
    g = create_group(f, groupname)
    for (name, set) in pairs(sets)
        g[name] = collect( set )
    end
end

function _save_cells!(f, grid::Grid{dim, C}) where {dim, C}
    g = create_group(f, celltypenames[C])
    g["nodeids"] = [ cell.nodes for cell in grid.cells ]
end

###################################################################################################
###################################################################################################
# TODO: implement, test, document

celltypenames = Dict([
    Line => "Line", QuadraticLine => "QuadraticLine",
    Triangle => "Triangle", QuadraticTriangle => "QuadraticTriangle", 
    Quadrilateral => "Quadrilateral", QuadraticQuadrilateral => "QuadraticQuadrilateral",
    Tetrahedron => "Tetrahedron",
    Hexahedron => "Hexahedron"
    ])

function _load_nodes(f)
    g = f["nodes"]
    coords = read(g, "coords")
    dim = length(coords[1])
    nodes = [ Node(Tensor{1, dim}(values(c))) for c in coords ]
    return nodes
end

function _load_sets(f, groupname::String, T::Type)
    g = f[groupname]
    if isempty(keys(g))
        pairs = Pair{String,Set{T}}[]
    elseif T <: Number
        pairs = Pair{String,Set{T}}[ key => Set{T}(read(g, key)) for key in keys(g) ]
    else
        pairs = Pair{String,Set{T}}[ key => Set{T}([ T(values(nt...)) for nt in read(g, key) ]) for key in keys(g) ]
    end
    return Dict{String,Set{T}}(pairs)
end

function _load_cells(f)
    cells = nothing
    for (C, cellname) in celltypenames
        if haskey(f, cellname)
            g = f[cellname]
            cells = collect( C(values(n)) for n in read(g, "nodeids") )
            return cells
        end
    end
    if isnothing(cells)
        throw(error("The file does not contain a group for one of the following cell types: $(values(celltypenames))"))
    end
    return cells
end

"""
    load(filepath::String)

Returns a grid saved usedin `@refpoint` from an HDF5 file saved at `filepath`.
"""
function load(filepath::String)
    f = h5open(filepath, "r")
    nodes = _load_nodes(f)
    facetsets = _load_sets(f, "facet sets", FacetIndex)
    cellsets  = _load_sets(f, "cell sets", Int)
    nodesets  = _load_sets(f, "node sets", Int)
    cells = _load_cells(f)
    close(f)
    return Grid(cells, nodes; facetsets=facetsets, cellsets=cellsets, nodesets=nodesets)
end