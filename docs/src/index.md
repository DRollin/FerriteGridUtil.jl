```@meta
CurrentModule = FerriteGridUtil
```

# FerriteGridUtil.jl

Welcome to the documentation for FerriteGridUtil!
This package adds functionality to the
finite element toolbox [Ferrite](https://github.com/Ferrite-FEM/Ferrite.jl) for working
with grids.

## How the documentation is organized

The documentation assumes that you are already familiar with the basic usage of Ferrite.
If not, you should first take a look at the [Ferrite documentation](https://ferrite-fem.github.io/Ferrite.jl/dev/).
Here, only the additional tools are explained.

After a basic introduction on this side, the document is organized as follows:

 - [**Tutorials**](tutorials/index.md) are documented examples.
 - [**Reference**](reference/index.md) contains the technical API reference of functions and
   methods (e.g. the documentation strings).

## Getting started

As a new user of FerriteGridUtil it is suggested to read the introduction on this side
and then start working with the tutorials before using FerriteGridUtil.

### Getting help

If you have questions about FerriteGridUtil it is suggested to use the `#ferrite-fem` channel on the
[Julia Slack](https://julialang.org/slack/), or the `#Ferrite.jl` stream on
[Zulip](https://julialang.zulipchat.com/).

## Important
The documentation below is just a template, FerriteGridUtil.jl is not a registered package!  

### Installation

To use FerriteGridUtil you first need to install Julia, see <https://julialang.org/> for details.
Installing FerriteGridUtil can then be done from the Pkg REPL; press `]` at the `julia>` promp to
enter `pkg>` mode:

```
pkg> add FerriteGridUtil
```

This will install FerriteGridUtil and all necessary dependencies. Press backspace to get back to the
`julia>` prompt. (See the [documentation for Pkg](https://pkgdocs.julialang.org/), Julia's
package manager, for more help regarding package installation and project management.)

Note, that you also need to install Ferrite, which can be done in the same way.

Finally, to load Ferrite and FerriteGridUtil, use

```julia
using Ferrite, FerriteGridUtil
```

You are now all set to start using FerriteGridUtil!

### Introduction to functionalities

