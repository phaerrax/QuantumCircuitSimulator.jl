module QuantumCircuitSimulator

using ITensors,
    LindbladVectorizedTensors,
    JSON,
    OffsetArrays,
    LinearAlgebra,
    RBNF,
    OpenQASM,
    DataStructures

include("utils.jl")

export gate, GateName, @GateName_str
include("gatename.jl")
include("gate.jl")

export PauliString, indices, operators, order, SPLNoiseModel, nqbits
include("paulistring.jl")

#include("qelib1_gates_2.0.jl")
#include("qelib1_gates_3.0.jl")
include("qiskit_gates.jl")
include("gates_vqbits.jl")

export compose, gates, gatelayers
include("qasm_itensors_parser.jl")

export noiselayer, inversenoiselayer
include("spl_noise_model.jl")

end
