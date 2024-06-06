using QuantumCircuitSimulator
using ITensors, LindbladVectorizedTensors
using OpenQASM
using Conda, PyCall, Pkg
using Test

Conda.add(["qiskit", "qiskit-aer"])
ENV["PYTHON"] = ""
Pkg.build("PyCall")

include("replicateqiskit.jl")

@testset "Qiskit compatibility" begin
    @test replicateqiskit(qiskitcircuit_noentanglement)
    @test replicateqiskit(qiskitcircuit_unitarygates)
    @test replicateqiskit(qiskitcircuit_registermapping)
    @test paulistringordering("XIYXIYZXII")
end

include("noise_layer_preserves_trace.jl")

@testset "Noise layer construction" begin
    @test noise_layer_preserves_trace(; N=20)
end

include("add_gate_from_file.jl")

@testset "Automatic gate definition" begin
    @test add_gate_from_file()
end

include("bell_state_from_openqasm.jl")

@testset "Build Bell state from complicated OpenQASM circuit" begin
    @test bell_state_from_openqasm(; atol=1e-12)
end
