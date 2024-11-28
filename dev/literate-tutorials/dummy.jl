# # [Dummy tutorial](@id tutorial-dummy)
#
# ## Introduction
#
# This is just a short dummy tutorial
#
# ## Commented Program
#
# What follows is a program spliced with comments.
#md # The full program, without comments, can be found in the next [section](@ref tutorial-dummy-plain-program).
#
# First we load all the packages we need for the computation.

using Ferrite, FerriteGridUtil

# Then, we do a computation.

1 + 1 == 2 ? true : false

#md # ## [Plain program](@id tutorial-dummy-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`dummy.jl`](dummy.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```