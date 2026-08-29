# Building and running circuits

A [`QuantumCircuit`](@ref) stores the qbit sites of a circuit together with its
sequence of gates, arranged in layers. Each layer is a list of gates that act on
disjoint qbits, so that all gates within a layer commute with each other and can
in principle be applied at the same time.

## Parsing an OpenQASM program

The most common way to build a circuit is to parse it directly from OpenQASM 2.0
code, either as a string or as an already-parsed `OpenQASM.Types.MainProgram`:

```jldoctest circuits
julia> using QuantumCircuitSimulator

julia> qasm = "OPENQASM 2.0; qreg q[4]; h q[0]; cx q[0], q[1]; cx q[1], q[2];";

julia> circ = QuantumCircuit(qasm);
```

Gates are assigned to layers following an approximation of Qiskit's "ASAP"
scheduling policy: a gate is placed in the earliest layer where all the qbits it
acts on are still free. Barrier instructions are taken into account: they
prevent gates on the barrier's qbits from being scheduled earlier than the
barrier.

The `depth` and `nqbits` functions return the number of layers and the number of
qbits in the circuit, respectively:

```jldoctest circuits
julia> depth(circ)  # number of layers
3

julia> nqbits(circ)  # number of qbits
4
```

A `QuantumCircuit` is iterable: iterating over it yields its layers (each one a
list of `ITensor` gates) in order.

### Operator picture

By default, the circuit is built to act on pure states, using the `"Qubit"` site
type from ITensors. Passing `operator_picture=true` builds the circuit using the
`"vQubit"` site type (from
[LindbladVectorizedTensors](https://github.com/phaerrax/LindbladVectorizedTensors.jl))
instead, appropriate for gates acting on vectorised density matrices or
observables (e.g. in the Heisenberg picture):

```jldoctest circuits
julia> QuantumCircuit(qasm; operator_picture=true);
```

### Custom gate definitions

OpenQASM programs may include custom `gate` declarations. These are parsed and
registered as new methods of the [`gate`](@ref) function as they are
encountered, so that later instructions in the same program can use them.

!!! warning "Overwriting existing gates"
    If a gate with the same name and the same signature already exists, it is
    overwritten by the new definition. This would trigger "method overwritten"
    warnings by the Julia compiler, that by default is suppressed; to show them,
    pass the `warn_on_gate_redefinition=true` to the `QuantumCircuit`
    constructor.

## From circuit to tensor network

Once a `QuantumCircuit` has been built, its layers can be converted to MPOs by
[`layers_mpo`](@ref):

```julia
julia> mpos = layers_mpo(circ);
```

The resulting list of MPOs can then be applied in sequence to an MPS (or to
another MPO) with ITensorMPS.jl's `apply` function to simulate the action of the
circuit.
