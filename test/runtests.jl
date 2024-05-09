using TEM
using ITensors
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
end
