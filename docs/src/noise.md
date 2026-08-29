# Sparse Pauli-Lindblad noise model

The [`SPLNoiseModel`](@ref) type represents a sparse Pauli-Lindblad (SPL) noise
model: a collection of one- and two-body Pauli-string Lindblad generators, each
with its own coefficient, acting on a specific set of qbits.  This is the noise
model described for example in [Filippov2023](@cite).

## Constructing a noise model

A model can be built directly from a dictionary mapping `PauliString`s (from
PauliStringTensors.jl) to their coefficients, together with the list of qbit
indices the strings refer to:

```jldoctest noise
julia> using PauliStringTensors, QuantumCircuitSimulator

julia> params = Dict("XI" => 0.01, "ZI" => 0.02, "XX" => 0.005);

julia> model = SPLNoiseModel(params, 1:2);

```

Only Pauli strings of order 1 or 2 (on adjacent positions) are allowed.
The second argument can be omitted, in which case the qbit range is assumed to
be `1:N`, with `N` the length of the Pauli strings. A model can also be loaded
from a JSON file mapping Pauli strings to coefficients, such as

```julia
model = SPLNoiseModel("noise_parameters.json")
```

## Noise layers

Given the ITensor site indices of a circuit (necessarily with the `"vQubit"`
site  type), [`noiselayer`](@ref) builds the MPO representation of one layer of
the noise channel described by an SPL model:

```jldoctest noise
julia> using ITensorMPS

julia> sites = siteinds("vQubit", 4);

julia> noise_mpo = noiselayer(sites, model);
```

An optional `maxbonddim` keyword argument truncates the resulting MPO to a
maximum bond dimension. The `MPO(sites, model::SPLNoiseModel)` constructor is an
alias for the same operation.

Finally, [`inversenoiselayer`](@ref) builds the MPO for the *inverse* of the
noise channel (useful for error mitigation), using the same interface as
`noiselayer`.
