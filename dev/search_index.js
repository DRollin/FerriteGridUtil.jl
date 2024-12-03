var documenterSearchIndex = {"docs":
[{"location":"devdocs/example/#Title","page":"Title","title":"Title","text":"","category":"section"},{"location":"devdocs/example/#Section-title","page":"Title","title":"Section title","text":"","category":"section"},{"location":"devdocs/example/","page":"Title","title":"Title","text":"","category":"page"},{"location":"reference/#Reference","page":"Reference overview","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference overview","title":"Reference overview","text":"Pages = [\n    \"property.md\",\n    \"manipulation.md\",\n    \"saveload.md\",\n    \"conversion.md\",\n]","category":"page"},{"location":"tutorials/#Tutorials","page":"Tutorials overview","title":"Tutorials","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials overview","title":"Tutorials overview","text":"On this page you find a collection of examples on how to use the provided functions for working with grids in Ferrite.If you are new to Ferrite, you should start there to learn the basics.","category":"page"},{"location":"tutorials/#[Tutorial-1:-Dummy](dummy.md)","page":"Tutorials overview","title":"Tutorial 1: Dummy","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials overview","title":"Tutorials overview","text":"This tutorial is just a dummy.","category":"page"},{"location":"tutorials/","page":"Tutorials overview","title":"Tutorials overview","text":"Keywords: dummy 1, dummy 2.","category":"page"},{"location":"devdocs/#Developer-documentation","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"devdocs/","page":"Developer documentation","title":"Developer documentation","text":"Here you can find some documentation of the internals of FerriteGridUtil.jl which are useful when developing the library.","category":"page"},{"location":"devdocs/","page":"Developer documentation","title":"Developer documentation","text":"Depth = 1\nPages = [    \n    \"example.md\",\n    ]","category":"page"},{"location":"reference/manipulation/#Manipulations","page":"Manipulations","title":"Manipulations","text":"","category":"section"},{"location":"reference/manipulation/","page":"Manipulations","title":"Manipulations","text":"Modules = [FerriteGridUtil]\nPages   = [\"manipulation.jl\"]","category":"page"},{"location":"reference/manipulation/#FerriteGridUtil.scale_relative!-Union{Tuple{dim}, Tuple{Grid{dim}, Real}} where dim","page":"Manipulations","title":"FerriteGridUtil.scale_relative!","text":"scale_relative!()\n\nScale the grid grid by a factor scalefactor relative to a point refpoint. The default refpoint is [0.0, 0.0, 0.0].\n\n\n\n\n\n","category":"method"},{"location":"reference/manipulation/#FerriteGridUtil.scale_relative-Union{Tuple{dim}, Tuple{Grid{dim}, Real}} where dim","page":"Manipulations","title":"FerriteGridUtil.scale_relative","text":"scale_relative(grid::Grid{dim}, scalefactor::Real; refpoint::Vec{dim}=zero(Vec{dim})) where {dim}\n\nReturn a grid scaled by scalefactor relative to a point refpoint. The default refpoint is [0.0, 0.0, 0.0].\n\n\n\n\n\n","category":"method"},{"location":"reference/manipulation/#FerriteGridUtil.shift_by!-Union{Tuple{dim}, Tuple{Grid{dim}, Vec{dim}}} where dim","page":"Manipulations","title":"FerriteGridUtil.shift_by!","text":"shift_by!(grid::Grid{dim}, shiftvalue::Vec{dim}) where {dim}\n\nShift the grid grid by a value shiftvalue.\n\n\n\n\n\n","category":"method"},{"location":"reference/manipulation/#FerriteGridUtil.shift_by-Union{Tuple{dim}, Tuple{Grid{dim}, Vec{dim}}} where dim","page":"Manipulations","title":"FerriteGridUtil.shift_by","text":"shift_by(grid::Grid{dim}, shiftvalue::Vec{dim}) where {dim}\n\nReturn a grid shifted by shiftvalue. \n\n\n\n\n\n","category":"method"},{"location":"reference/conversion/#Converting-to-other-mesh-types","page":"Converting to other mesh types","title":"Converting to other mesh types","text":"","category":"section"},{"location":"reference/conversion/","page":"Converting to other mesh types","title":"Converting to other mesh types","text":"Modules = [FerriteGridUtil]\nPages   = [\"conversion.jl\"]","category":"page"},{"location":"reference/conversion/#FerriteGridUtil.convert_to_makie_mesh-Union{Tuple{dim}, Tuple{Grid{dim}, Union{AbstractSet{Int64}, String}}} where dim","page":"Converting to other mesh types","title":"FerriteGridUtil.convert_to_makie_mesh","text":"convert_to_makie_mesh()\n\nReturn a mesh that can be used for plotting with Makie.jl.\n\n\n\n\n\n","category":"method"},{"location":"reference/saveload/#Saving-and-loading-to/from-a-file","page":"Saving and loading to/from a file","title":"Saving and loading to/from a file","text":"","category":"section"},{"location":"reference/saveload/","page":"Saving and loading to/from a file","title":"Saving and loading to/from a file","text":"Modules = [FerriteGridUtil]\nPages   = [\"saveload.jl\"]","category":"page"},{"location":"reference/saveload/#FerriteGridUtil.load-Tuple{String}","page":"Saving and loading to/from a file","title":"FerriteGridUtil.load","text":"load(filepath::String)\n\nReturns a grid saved usedin @refpoint from an HDF5 file saved at filepath.\n\n\n\n\n\n","category":"method"},{"location":"reference/saveload/#FerriteGridUtil.save-Union{Tuple{C}, Tuple{dim}, Tuple{Grid{dim, C}, String}} where {dim, C}","page":"Saving and loading to/from a file","title":"FerriteGridUtil.save","text":"save(grid::Grid{dim, C}, filepath::String) where {dim, C}\n\nSaves a grid grid in an HDF5 file at filepath.\n\n\n\n\n\n","category":"method"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"EditURL = \"../literate-tutorials/dummy.jl\"","category":"page"},{"location":"tutorials/dummy/#tutorial-dummy","page":"Dummy tutorial","title":"Dummy tutorial","text":"","category":"section"},{"location":"tutorials/dummy/#Introduction","page":"Dummy tutorial","title":"Introduction","text":"","category":"section"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"This is just a short dummy tutorial","category":"page"},{"location":"tutorials/dummy/#Commented-Program","page":"Dummy tutorial","title":"Commented Program","text":"","category":"section"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"What follows is a program spliced with comments. The full program, without comments, can be found in the next section.","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"First we load all the packages we need for the computation.","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"using Ferrite, FerriteGridUtil","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"Then, we do a computation.","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"1 + 1 == 2 ? true : false","category":"page"},{"location":"tutorials/dummy/#tutorial-dummy-plain-program","page":"Dummy tutorial","title":"Plain program","text":"","category":"section"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"Here follows a version of the program without any comments. The file is also available here: dummy.jl.","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"using Ferrite, FerriteGridUtil\n\n1 + 1 == 2 ? true : false","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"","category":"page"},{"location":"tutorials/dummy/","page":"Dummy tutorial","title":"Dummy tutorial","text":"This page was generated using Literate.jl.","category":"page"},{"location":"reference/property/#Computing-properties","page":"Computing properties","title":"Computing properties","text":"","category":"section"},{"location":"reference/property/","page":"Computing properties","title":"Computing properties","text":"Modules = [FerriteGridUtil]\nPages   = [\"property.jl\"]","category":"page"},{"location":"reference/property/#FerriteGridUtil.get_coordinate_limits-Union{Tuple{Grid{dim}}, Tuple{dim}} where dim","page":"Computing properties","title":"FerriteGridUtil.get_coordinate_limits","text":"get_coordinate_limits(grid::Grid{dim}) where {dim}\n\nReturn a tuple of the minimum and maximum values for each grid direction.\n\n\n\n\n\n","category":"method"},{"location":"reference/property/#FerriteGridUtil.get_dofs_from_coordinate-Union{Tuple{dim}, Tuple{DofHandler{dim, G} where G<:Ferrite.AbstractGrid{dim}, Vec{dim}, Symbol}} where dim","page":"Computing properties","title":"FerriteGridUtil.get_dofs_from_coordinate","text":"get_dofs_from_coordindate(dh::DofHandler{dim}, x::Vec{dim}, fieldname::Symbol; radius::Real=1e-3) where {dim}\n\nReturn the degrees of freedom corresponding to fieldname for a node at coordinate x.  The node must be within a neighbourhood of radius radius. The default radius is 1e-3.\n\n\n\n\n\n","category":"method"},{"location":"reference/property/#FerriteGridUtil.get_interface_between_sets-Union{Tuple{dim}, Tuple{Grid{dim}, AbstractSet{Int64}, AbstractSet{Int64}}} where dim","page":"Computing properties","title":"FerriteGridUtil.get_interface_between_sets","text":"get_interface_between_sets(grid::Grid{dim}, set¹::AbstractSet{Int}, set²::AbstractSet{Int}) where {dim}\n\nReturn an OrderedSet{FacetIndex} of the facets building an interface between the two sets set¹ and set².\n\n\n\n\n\n","category":"method"},{"location":"reference/property/#FerriteGridUtil.get_moment-Union{Tuple{dim}, Tuple{Grid{dim}, Int64}} where dim","page":"Computing properties","title":"FerriteGridUtil.get_moment","text":"get_moment()\n\nReturn the moment of a grid or a cellset. The order ∈ {0, 1, 2} correspondes to the  volume, centre of gravity and second moment of inertia respectively. The moment is calculated using the default QuadratureRule for the cells of the grid.\n\n\n\n\n\n","category":"method"},{"location":"reference/property/#FerriteGridUtil.get_node_from_coordinate-Union{Tuple{dim}, Tuple{DofHandler{dim, G} where G<:Ferrite.AbstractGrid{dim}, Vec{dim}}} where dim","page":"Computing properties","title":"FerriteGridUtil.get_node_from_coordinate","text":"get_node_from_coordinate(dh::DofHandler{dim}, x::Vec{dim}; radius::Real=1e-3) where {dim}\n\nSearch grid for the first node within a neighbourhood of radius radius around x. The node id, a cell id to which the node belongs and a bool to signify a succesful search are returned.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FerriteGridUtil","category":"page"},{"location":"#FerriteGridUtil.jl","page":"Home","title":"FerriteGridUtil.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for FerriteGridUtil! This package adds functionality to the finite element toolbox Ferrite for working with grids.","category":"page"},{"location":"#How-the-documentation-is-organized","page":"Home","title":"How the documentation is organized","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The documentation assumes that you are already familiar with the basic usage of Ferrite. If not, you should first take a look at the Ferrite documentation. Here, only the additional tools are explained.","category":"page"},{"location":"","page":"Home","title":"Home","text":"After a basic introduction on this side, the document is organized as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Tutorials are documented examples.\nReference contains the technical API reference of functions and methods (e.g. the documentation strings).","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As a new user of FerriteGridUtil it is suggested to read the introduction on this side and then start working with the tutorials before using FerriteGridUtil.","category":"page"},{"location":"#Getting-help","page":"Home","title":"Getting help","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you have questions about FerriteGridUtil it is suggested to use the #ferrite-fem channel on the Julia Slack, or the #Ferrite.jl stream on Zulip.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To use FerriteGridUtil you first need to install Julia, see https://julialang.org/ for details. Installing FerriteGridUtil can then be done from the Pkg REPL; press ] at the julia> promp to enter pkg> mode:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add FerriteGridUtil","category":"page"},{"location":"","page":"Home","title":"Home","text":"This will install FerriteGridUtil and all necessary dependencies. Press backspace to get back to the julia> prompt. (See the documentation for Pkg, Julia's package manager, for more help regarding package installation and project management.)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note, that you also need to install Ferrite, which can be done in the same way.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, to load Ferrite and FerriteGridUtil, use","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Ferrite, FerriteGridUtil","category":"page"},{"location":"","page":"Home","title":"Home","text":"You are now all set to start using FerriteGridUtil!","category":"page"},{"location":"#Introduction-to-functionalities","page":"Home","title":"Introduction to functionalities","text":"","category":"section"}]
}
