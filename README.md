# QuantumCircuitSimulator.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://phaerrax.github.io/QuantumCircuitSimulator.jl/dev/)

QuantumCircuitSimulator.jl loads quantum circuits from
[OpenQASM](https://openqasm.com/) code and simulates them as tensor networks,
using the [ITensor](https://itensor.org/) library. Circuits are represented as
sequences of layers of gates, which can then be turned into matrix product
operators (MPOs) and applied to matrix product states (MPSs), enabling classical
simulation of circuits that would otherwise be too large to represent as dense
state vectors.

The package also provides tools to describe and apply a sparse Pauli-Lindblad
noise model, and to inspect Pauli-string decompositions of matrix product states
(useful, for example, when representing vectorized density matrices or
observables).

Documentation and installation instructions are available
[here](https://phaerrax.github.io/QuantumCircuitSimulator.jl/dev/).
