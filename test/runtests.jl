using TEM
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
