module FerriteGridUtil
using Ferrite
using OrderedCollections: OrderedSet
using Makie: Makie, GeometryBasics
using HDF5: h5open, create_group, read, close

# Internal helper functions
include("helper.jl")

# Functions to compute properties of a grid
include("property.jl")
export get_moment, get_coordinate_limits, get_interface_between_sets

# Functions to manipulate a grid
include("manipulation.jl")
export scale_relative, scale_relative!, shift_by, shift_by!

# Functions to save and load a grid
include("saveload.jl")
export save, load

# Functions to convert a grid to a mesh-type of a different package
include("conversion/conversion.jl")
export convert_to_makie_mesh, disconnect_field

end
