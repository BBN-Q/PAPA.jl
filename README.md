![PAPA](docs/images/PAPA_BEAR_small.png)

**PAPA BEAR**: **PA**irwise **P**erturbative **A**nsatz for **B**etter **E**stimation **A**nd **R**untime written in Julia

[![CI](https://github.com/BBN-Q/PAPA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/BBN-Q/PAPA.jl/actions/workflows/CI.yml)

PAPA.jl is a Julia package used for bootstrapping the process matrix of a large quantum system from the descriptions of many smaller quantum systems.  In particular, PAPA.jl takes the processes matrices for all pairs of qubits in a N-qubit system and bootstraps to a process matrix for the full system. It does so in a way that constrains the global process matrix based on __all__ the pairwise reconstructions with a non-linear least-squares method.

## References
[arXiv](https://arxiv.org/abs/1902.10821)

[nature comm.](https://www.nature.com/articles/s41467-020-14873-1)

## Installation
This branch can now be locally installed as a package with the command:
```julia
(v1.0) pkg> add https://github.com/BBN-Q/PAPA.jl
```

To run the test code do:
```julia
(v1.0) pkg> test PAPA
```

## Usage
`papa_reconstruction()` is the primary function for PAPA reconstructions.

## License
Apache License v2.0

## Contributors
Luke Govia, Guilhem Ribeill, Matthew Ware, and Hari Krovi

## Acknowledgements
This effort was partially funded by ARO under contract W911NF-14-C-0048.

## Copyright
Raytheon BBN Technologies
