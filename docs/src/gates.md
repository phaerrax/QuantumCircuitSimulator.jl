# Gates

Quantum gates are represented as plain `ITensor` objects, built through the
[`gate`](@ref) function. This mirrors ITensors.jl's own `op` mechanism, but is
specialized to the needs of this package, in order to have a structured
way of building them from the OpenQASM syntax.

## Building a gate

The main entry point is the following method:

```julia
gate(name::AbstractString, s::Index...; cargs)
```

which returns the ITensor for the gate `name` acting on the given site indices
`s`. If the gate takes numerical parameters (angles, for example), they are
passed in the `cargs` keyword argument, in the same order used by the
OpenQASM/Qiskit definition of the gate:

```jldoctest gates; setup = :(using QuantumCircuitSimulator)
julia> using ITensorMPS

julia> s = siteinds("Qubit", 2);

julia> gate("x", s[1]);

julia> gate("u3", s[1]; cargs=(0, pi / 4, 2));

julia> gate("cx", s[1], s[2]);

```

A vector of sites together with integer positions can be used instead of the site
indices directly:

```jldoctest gates
julia> s = siteinds("Qubit", 4);

julia> gate("x", s, 2);

julia> gate("cnot", s, 1, 3);

```

Gates are resolved by dispatching on the common `SiteType` of the given indices
and on a [`GateName`](@ref), a parametrized type built from the gate name
string, like ITensors does with `OpName` or `StateName`. 
The gate format appropriate to pure states or vectorised density matrices is
automatically chosen depending on the site type:

* if the indices have the `"Qubit"` site type, then the standard version will be
  returned, where a gate \\(g\\) acts on a pure state \\(\psi\\) as \\(g\psi\\);
* if the indices have the `"vQubit"` site type, then the vectorised version will
  be returned, where a gate \\(g\\) acts on a mixed state \\(\rho\\) as
  \\(g\rho\adj{g}\\);

With this mechanism, new gates can be defined simply by adding new `gate`
methods (see below). New gates must be define only through the
`::SiteType"Qubit"` method: `::SiteType"vQubit"` gates are automatically
constructed from the corresponding pure-state versions.

## Defining new gates

A gate for a specific site type is defined by adding a method to `gate` that
dispatches on `GateName"..."`, with `SiteType"Qubit"`: for example

```julia
function QuantumCircuitSimulator.gate(::GateName"myz", ::SiteType"Qubit", s::Index)
    return gate("z", s)
end
```

Alternatively, a method returning a plain Julia matrix (rather than an
`ITensor`) can be defined by dispatching on `GateName` and `SiteType` alone,
without any index arguments); `gate` will then build the corresponding operator
automatically:

```julia
function QuantumCircuitSimulator.gate(::GateName"myx", ::SiteType"Qubit")
    return [0 1; 1 0]
end
```

The package already defines a standard library of gates (from the
"Quantum Experience" OpenQASM 2.0 header, plus a few extras) for both the
`"Qubit"` site type (pure-state simulation) and the `"vQubit"` site type
(vectorised states/operator simulation).
Custom gates declared in an OpenQASM program are registered the same way,
automatically, while the program is being parsed (see [Building and running
circuits](@ref)), so there is no need to define them manually before running the
simulation.

## Building a gate from a matrix directly

If a gate is only needed once, it does not have to be registered as a named
method: an ITensor can be built directly from a matrix and a set of site indices
with

```julia
gate(M::AbstractArray, s::Index...)
```

For example:

```jldoctest gates
julia> s = siteind("Qubit");

julia> g = gate([1/2 0; 0 -1/2], s);
```
