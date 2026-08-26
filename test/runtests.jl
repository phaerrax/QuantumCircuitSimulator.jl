using QuantumCircuitSimulator
using ITensors, ITensorMPS, LindbladVectorizedTensors, PauliStringTensors
using OpenQASM
using Conda, PyCall, Pkg
using Test

Conda.add(["qiskit", "qiskit-aer"])
ENV["PYTHON"] = ""
Pkg.build("PyCall")

include("utils.jl")

include("replicateqiskit.jl")

@testset "Qiskit compatibility" begin
    @test replicateqiskit(qiskitcircuit_noentanglement)
    @test replicateqiskit(qiskitcircuit_unitarygates)
    @test replicateqiskit(qiskitcircuit_registermapping)
    @test paulistringordering("XIYXIYZXII")
end

include("noise_layer_preserves_trace.jl")
include("spl_noise_truncation.jl")

@testset "Noise layer construction" begin
    @test noise_layer_preserves_trace(; N=20)
    @test spl_noise_truncation()
end

include("add_gate_from_file.jl")

@testset "Automatic gate definition" begin
    @test add_gate_from_file()
end

include("bellstate.jl")

@testset "Build Bell state from complicated OpenQASM circuit" begin
    @test bellstate_openqasm()
end

include("gates.jl")

@testset "Gate adjoints and Qubit/vQubit equivalence" begin
    test_Ugates()
    test_xyz_gates()
    test_hspt_gates()
    test_rotation_gates()
    test_other_gates()
    test_controlled_gates()
    test_gate_adjoints()
    test_gates_vqubit()
    test_ecr_vqubit()
end

include("pauli_sampling.jl")

@testset "Sampling of Pauli strings from MPS" begin
    @test pauli_sampling()
end

include("circuit_parsing.jl")

@testset "Circuit parsing" begin
    @test test_parse()
end
