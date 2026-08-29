# QuantumCircuitSimulator.jl

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

## Installation

### From a registry

The package is not registered in the General registry, but in my
[TensorNetworkSimulations](https://github.com/phaerrax/TensorNetworkSimulations)
registry. Adding this registry first, with

```julia
using Pkg
Pkg.Registry.add("https://github.com/phaerrax/TensorNetworkSimulations.git")
```

(this must be done just once per Julia installation) then install the package
as you would with a normal one:

```julia
using Pkg
Pkg.add("QuantumCircuitSimulator")
```

### From GitHub

Alternatively, straight installation from GitHub is also possible:

```julia
using Pkg
Pkg.add("https://github.com/phaerrax/QuantumCircuitSimulator.jl")
```

## Overview

* [Building and running circuits](@ref) explains how to load an OpenQASM program
  into a [`QuantumCircuit`](@ref) and turn it into a sequence of MPOs.
* [Gates](@ref) describes the `gate` interface used to build and extend the
  library of quantum gates.
* [Sparse Pauli-Lindblad noise model](@ref) covers the [`SPLNoiseModel`](@ref)
  type and the associated noise (and inverse-noise) layers.
* [Library reference](@ref) lists the available public methods.

## Quick example

```jldoctest
julia> using QuantumCircuitSimulator

julia> qasm = "OPENQASM 2.0; qreg q[2]; h q[0]; cx q[0], q[1];";

julia> circ = QuantumCircuit(qasm);

julia> mpos = layers_mpo(circ);

```

## Bibliography

```@bibliography
```
