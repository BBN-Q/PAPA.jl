![PAPA](docs/images/PAPA_BEAR_small.png)

**PA**irwise **P**erturbative **A**nsatz for **B**etter **E**stimation **A**nd **R**untime written in Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://matthewware.github.io/PAPA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://matthewware.github.io/PAPA.jl/dev)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://matthewware.gitlab.io/PAPA.jl/dev)
[![Build Status](https://travis-ci.com/matthewware/PAPA.jl.svg?branch=master)](https://travis-ci.com/matthewware/PAPA.jl)
[![Build Status](https://gitlab.com/matthewware/PAPA.jl/badges/master/build.svg)](https://gitlab.com/matthewware/PAPA.jl/pipelines)
[![Coverage](https://gitlab.com/matthewware/PAPA.jl/badges/master/coverage.svg)](https://gitlab.com/matthewware/PAPA.jl/commits/master)

PAPA.jl is a Julia package used for bootstrapping the process matrix of a large quantum system from the descriptions of many smaller quantum systems.  In particular, PAPA.jl takes a process matrix of an N dimensional system and bootstrapps to a process matrix for an M > N dimensional system but does so in a way that constrains the global process matrix based on __all__ the local reconstructions with a non-linear least-squares method.     

## Installation
This branch can now be locally installed as a package with the command:
```julia
(v1.0) pkg> add <path to this repo>
```

To run the test code do:
```julia
(v1.0) pkg> test PAPA
```

## Usage
`papa_reconstruction()` is the primary function for PAPA reconstructions.
