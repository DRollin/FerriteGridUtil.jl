_convert_cell_to_makie(cell::Union{Line,QuadraticLine}) = tuple(Makie.GeometryBasics.NgonFace{2,Int}(Ferrite.vertices(cell)))
_get_makie_type_data(::Union{Type{Line},Type{QuadraticLine}}) = 1, Makie.GeometryBasics.NgonFace{2,Int}

_convert_cell_to_makie(cell::Union{Triangle,QuadraticTriangle}) = tuple(Makie.GeometryBasics.NgonFace{3,Int}(Ferrite.vertices(cell)))
_get_makie_type_data(::Union{Type{Triangle},Type{QuadraticTriangle}}) = 1, Makie.GeometryBasics.NgonFace{3,Int}

_convert_cell_to_makie(cell::Union{Quadrilateral,QuadraticQuadrilateral}) = tuple(Makie.GeometryBasics.NgonFace{4,Int}(Ferrite.vertices(cell)))
_get_makie_type_data(::Union{Type{Quadrilateral},Type{QuadraticQuadrilateral}}) = 1, Makie.GeometryBasics.NgonFace{4,Int}

_convert_cell_to_makie(cell::Union{Tetrahedron,QuadraticTetrahedron}) = Tuple(Makie.GeometryBasics.NgonFace{3,Int}(facet) for facet in Ferrite.facets(cell))
_get_makie_type_data(::Union{Type{Tetrahedron},Type{QuadraticTetrahedron}}) = 4, Makie.GeometryBasics.NgonFace{3,Int}

_convert_cell_to_makie(cell::Union{Hexahedron,QuadraticHexahedron}) = Tuple(Makie.GeometryBasics.NgonFace{4,Int}(facet) for facet in Ferrite.facets(cell))
_get_makie_type_data(::Union{Type{Hexahedron},Type{QuadraticHexahedron}}) = 6, Makie.GeometryBasics.NgonFace{4,Int}

_convert_vec_to_makie(v::Vec{dim}) where {dim} = Makie.GeometryBasics.Point{dim, Float64}(v...)