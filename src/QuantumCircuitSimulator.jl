module QuantumCircuitSimulator

using DataStructures,
    ITensors,
    ITensorMPS,
    JSON,
    LindbladVectorizedTensors,
    LinearAlgebra,
    Memoize,
    OffsetArrays,
    OpenQASM,
    ProgressMeter,
    PauliStringTensors,
    RBNF

include("utils.jl")

export gate, GateName, @GateName_str
include("gatename.jl")
include("gate.jl")

export samplepaulistrings, relevantpaulistrings
include("paulistring_sampling.jl")

#include("qelib1_gates_2.0.jl")
#include("qelib1_gates_3.0.jl")
include("qiskit_gates.jl")
include("gates_vqbits.jl")

export compose, gates
include("qasm_itensors_parser.jl")

export QuantumCircuit, depth, quantumcircuit, parsecircuit, layers_mpo
include("quantum_circuit.jl")

export qbitsites, noiselayer, inversenoiselayer, crop, SPLNoiseModel, nqbits
include("spl_noise_model.jl")

export gatelayers
include("deprecated.jl")

end
