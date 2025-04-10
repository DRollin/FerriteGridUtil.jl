# TODO: test, document
"""
    scale_relative(grid::Grid, scalefactor::Real; refpoint::Vec)

Return a grid scaled by `scalefactor` relative to a point `refpoint`. The default `refpoint` is the `[0.0, 0.0, 0.0]`.
"""
function scale_relative(grid::Grid{dim}, scalefactor::Real; refpoint::Vec{dim}=zero(Vec{dim})) where {dim}
    return scale_relative!(deepcopy(grid), scalefactor; refpoint=refpoint)
end

"""
    scale_relative!(grid::Grid, scalefactor::Real; refpoint::Vec)

Scale the grid `grid` by a factor `scalefactor` relative to a point `refpoint`. The default `refpoint` is `[0.0, 0.0, 0.0]`.
"""
function scale_relative!(grid::Grid{dim}, scalefactor::Real; refpoint::Vec{dim}=zero(Vec{dim})) where {dim}
    for (i, n) in pairs(grid.nodes)
        grid.nodes[i] = Node(refpoint + scalefactor*(n.x - refpoint))
    end
    return grid
end

###################################################################################################
###################################################################################################
#TODO: test, document
"""
    shift_by(grid::Grid, shiftvalue::Vec)

Return a grid shifted by `shiftvalue`. 
"""
function shift_by(grid::Grid{dim}, shiftvalue::Vec{dim}) where {dim}
    return shift_by!(deepcopy(grid), shiftvalue)
end

"""
    shift_by!(grid::Grid, shiftvalue::Vec)

Shift the grid `grid` by a value `shiftvalue`.
"""
function shift_by!(grid::Grid{dim}, shiftvalue::Vec{dim}) where {dim}
    for (i, n) in pairs(grid.nodes)
        grid.nodes[i] = Node(n.x + shiftvalue)
    end
    return grid
end

