using LinearAlgebra

function isunitary(g::ITensor)
    site_inds = inds(g; plev=0)
    id = ITensors.op(I, site_inds...)
    return apply(gateadjoint(g), g) ≈ id
end

function rand_qubit_state(s::Index)
    return apply(ITensors.op("RandomUnitary", s), ITensors.state("0", s))
end

gateadjoint(g) = swapprime(dag(g), 0 => 1)

# Test that we defined the gates coming from the Qiskit library correctly.
# The following properties come from Qiskit's documentation (plus something more), so it's
# crucial that the gate satisfy them.
# This way we also make sure that each gate gets called and tested at least once.

function test_Ugates()
    q = siteind("Qubit")
    θ, ϕ, λ = rand(3)

    @test isunitary(gate("u1", q; cargs=(θ,)))
    @test isunitary(gate("u2", q; cargs=(θ, ϕ)))
    @test isunitary(gate("u3", q; cargs=(θ, ϕ, λ)))

    @test gate("u1", q; cargs=(-θ,)) ≈ gateadjoint(gate("u1", q; cargs=(θ,)))
    @test gate("u2", q; cargs=(pi-ϕ, pi-θ)) ≈ gateadjoint(gate("u2", q; cargs=(θ, ϕ)))
    @test gate("u3", q; cargs=(-θ, -λ, -ϕ)) ≈ gateadjoint(gate("u3", q; cargs=(θ, ϕ, λ)))

    @test gate("u3", q; cargs=(θ, ϕ, λ)) ≈ gate("u", q; cargs=(θ, ϕ, λ))
    @test gate("u3", q; cargs=(θ, ϕ, λ)) ≈
        cis((ϕ+λ)/2) * apply(
        gate("rz", q; cargs=(ϕ,)),
        apply(gate("ry", q; cargs=(θ,)), gate("rz", q; cargs=(λ,))),
    )
    @test gate("u3", q; cargs=(θ, ϕ, λ)) ≈ apply(
        gate("p", q; cargs=(ϕ,)), apply(gate("ry", q; cargs=(θ,)), gate("p", q; cargs=(λ,)))
    )
    @test gate("u3", q; cargs=(θ, -pi/2, pi/2)) ≈ gate("rx", q; cargs=(θ,))
    @test gate("u3", q; cargs=(θ, 0, 0)) ≈ gate("ry", q; cargs=(θ,))

    @test gate("u2", q; cargs=(θ, ϕ)) ≈ gate("u", q; cargs=(pi/2, θ, ϕ))
    @test gate("u2", q; cargs=(0, pi)) ≈ gate("h", q)
    @test gate("u2", q; cargs=(0, 0)) ≈ gate("ry", q; cargs=(pi/2,))
    @test gate("u2", q; cargs=(-pi/2, pi/2)) ≈ gate("rx", q; cargs=(pi/2,))
    @test gate("u2", q; cargs=(θ, ϕ)) ≈
        cis((θ+ϕ)/2) * apply(
        gate("rz", q; cargs=(θ,)),
        apply(gate("ry", q; cargs=(pi/2,)), gate("rz", q; cargs=(ϕ,))),
    )
    @test gate("u2", q; cargs=(θ, ϕ)) ≈
        cispi(-1/4) * apply(
        gate("p", q; cargs=(θ+pi/2,)), apply(gate("sx", q), gate("p", q; cargs=(ϕ-pi/2,)))
    )

    @test gate("u1", q; cargs=(θ,)) ≈ gate("u", q; cargs=(0, 0, θ))
    @test gate("u1", q; cargs=(pi,)) ≈ gate("z", q)
    @test gate("u1", q; cargs=(pi/2,)) ≈ gate("s", q)
    @test gate("u1", q; cargs=(pi/4,)) ≈ gate("t", q)
    @test gate("u1", q; cargs=(θ,)) ≈ cis(θ/2) * gate("rz", q; cargs=(θ,))
end

function test_xyz_gates()
    q = siteind("Qubit")

    @test gate("x", q) ≈ im * gate("rx", q; cargs=(pi,))
    @test apply(gate("x", q), state("0", q)) ≈ state("1", q)
    @test apply(gate("x", q), state("1", q)) ≈ state("0", q)

    @test gate("y", q) ≈ im * gate("ry", q; cargs=(pi,))
    @test apply(gate("y", q), state("0", q)) ≈ im * state("1", q)
    @test apply(gate("y", q), state("1", q)) ≈ -im * state("0", q)

    @test gate("z", q) ≈ im * gate("rz", q; cargs=(pi,))
    @test apply(gate("z", q), state("0", q)) ≈ state("0", q)
    @test apply(gate("z", q), state("1", q)) ≈ -state("1", q)
end

function test_hspt_gates()
    q = siteind("Qubit")
    θ = rand()

    @test apply(gate("h", q), state("0", q)) ≈ state("+", q)
    @test apply(gate("h", q), state("1", q)) ≈ state("-", q)
    @test apply(gate("h", q), state("+", q)) ≈ state("0", q)
    @test apply(gate("h", q), state("-", q)) ≈ state("1", q)

    @test apply(gate("s", q), gate("s", q)) ≈ gate("z", q)
    # S is equal to a pi/2 rotation around Z up to a pi/4 phase.
    @test gate("s", q) ≈ cispi(1/4)*gate("rz", q; cargs=(pi/2,))
    @test gate("sdg", q) ≈ cispi(-1/4)*gate("rz", q; cargs=(-pi/2,))
    @test gateadjoint(gate("s", q)) ≈ gate("sdg", q)

    @test gate("p", q; cargs=(θ,)) ≈ cis(θ/2) * gate("rz", q; cargs=(θ,))
    @test gate("p", q; cargs=(pi,)) ≈ gate("z", q)
    @test gate("p", q; cargs=(pi/2,)) ≈ gate("s", q)
    @test gate("p", q; cargs=(pi/4,)) ≈ gate("t", q)

    @test apply(gate("t", q), apply(gate("t", q), apply(gate("t", q), gate("t", q)))) ≈
        gate("z", q)
    # T is equal to a pi/4 rotation around Z up to a pi/8 phase.
    @test gate("t", q) ≈ cispi(1/8) * gate("rz", q; cargs=(pi/4,))
    @test gateadjoint(gate("t", q)) ≈ gate("tdg", q)
end

function test_rotation_gates()
    q = siteind("Qubit")
    θ = rand()

    @test gate("rx", q; cargs=(θ,)) ≈ exp(-im*θ/2 * gate("x", q))
    @test gate("ry", q; cargs=(θ,)) ≈ exp(-im*θ/2 * gate("y", q))
    @test gate("rz", q; cargs=(θ,)) ≈ exp(-im*θ/2 * gate("z", q))

    q1, q2 = siteinds("Qubit", 2)
    @test gate("rxx", q1, q2; cargs=(θ,)) ≈ exp(-im*θ/2 * gate("x", q1) * gate("x", q2))
    @test gate("rxx", q1, q2; cargs=(pi,)) ≈ -im * gate("x", q1) * gate("x", q2)
    @test gate("rxx", q1, q2; cargs=(θ,)) ≈ apply(
        gate("h", q1) * gate("h", q2),
        apply(
            gate("cx", q1, q2),
            apply(
                gate("rz", q2; cargs=(θ,)),
                apply(gate("cx", q1, q2), gate("h", q1) * gate("h", q2)),
            ),
        ),
    )

    @test gate("ryy", q1, q2; cargs=(θ,)) ≈ exp(-im*θ/2 * gate("y", q1) * gate("y", q2))
    @test gate("ryy", q1, q2; cargs=(θ,)) ≈ apply(
        gate("sxdg", q1) * gate("sxdg", q2),
        apply(
            gate("cx", q1, q2),
            apply(
                gate("rz", q2; cargs=(θ,)),
                apply(gate("cx", q1, q2), gate("sx", q1) * gate("sx", q2)),
            ),
        ),
    )

    @test gate("rzz", q1, q2; cargs=(θ,)) ≈ exp(-im*θ/2 * gate("z", q1) * gate("z", q2))
    @test gate("rzz", q1, q2; cargs=(pi,)) ≈ -im * gate("z", q1) * gate("z", q2)
    @test gate("rzz", q1, q2; cargs=(θ,)) ≈
        apply(gate("cx", q1, q2), apply(gate("rz", q2; cargs=(θ,)), gate("cx", q1, q2)))
end

function test_other_gates()
    q = siteinds("Qubit", 5)

    @test gate("sx", q[1]) ≈ cispi(1/4) * gate("rx", q[1]; cargs=(pi/2,))
    # Qiskit gets this wrong, in its manual it says 1/4 instead or -1/4.
    @test gate("sxdg", q[1]) ≈ cispi(-1/4) * gate("rx", q[1]; cargs=(-pi/2,))

    m1 = matrix(ITensors.op("RandomUnitary", siteind("Qubit")))
    m2 = matrix(ITensors.op("RandomUnitary", siteind("Qubit")))
    v00 = ITensors.state("0", q[1]) * ITensors.state("0", q[2])
    v12 = apply(ITensors.op(m1, q[1]) * ITensors.op(m2, q[2]), v00)
    v21 = apply(ITensors.op(m2, q[1]) * ITensors.op(m1, q[2]), v00)
    @test apply(gate("swap", q[1], q[2]), v12) ≈ v21

    # Note that `apply([A,B,C], X)` applies the elements in the list in the given order,
    # i.e. the result is `C(B(A(X)))`.
    @test gate("rccx", q[1], q[2], q[3]) ≈ apply(
        [
            gate("h", q[3]),
            gate("t", q[3]),
            gate("cx", q[2], q[3]),
            gate("tdg", q[3]),
            gate("cx", q[1], q[3]),
            gate("t", q[3]),
            gate("cx", q[2], q[3]),
            gate("tdg", q[3]),
            gate("h", q[3]),
        ],
        ITensors.op("Id", q[1], q[2], q[3]),
    )

    @test gate("rc3x", q[1], q[2], q[3], q[4]) ≈ apply(
        [
            gate("h", q[4]),
            gate("t", q[4]),
            gate("cx", q[3], q[4]),
            gate("tdg", q[4]),
            gate("h", q[4]),
            gate("cx", q[1], q[4]),
            gate("t", q[4]),
            gate("cx", q[2], q[4]),
            gate("tdg", q[4]),
            gate("cx", q[1], q[4]),
            gate("t", q[4]),
            gate("cx", q[2], q[4]),
            gate("tdg", q[4]),
            gate("h", q[4]),
            gate("t", q[4]),
            gate("cx", q[3], q[4]),
            gate("tdg", q[4]),
            gate("h", q[4]),
        ],
        ITensors.op("Id", q[1], q[2], q[3], q[4])
    )
end

# Create a controlled version of a gate, given a set of control qbits.
# TODO Can we promote this to a proper function in the gate library? It might come in handy
# to define all sorts of controlled gates in an easy way.
function control(gate::ITensor, control_qbit)
    @assert !in(control_qbit, inds(gate))
    return ITensors.op("Proj1", control_qbit) * gate +
           ITensors.op("Proj0", control_qbit) * ITensors.op(I, inds(gate; plev=0)...)
end

# Multi-site version: gate control is recursive.
function control(gate::ITensor, control_qbit1, control_qbits...)
    return control(control(gate, control_qbit1), control_qbits...)
end

function test_controlled_gates()
    q = siteinds("Qubit", 5)
    s0r = state("0", q[1]) * rand_qubit_state(q[2])
    s1r = state("1", q[1]) * rand_qubit_state(q[2])
    @test apply(gate("cx", q[1], q[2]), s0r) ≈ s0r
    @test apply(gate("cx", q[1], q[2]), s1r) ≈ apply(gate("x", q[2]), s1r)

    θ, ϕ, λ = rand(3)
    @test gate("cx", q[1], q[2]) ≈ control(gate("x", q[2]), q[1])
    @test gate("cnot", q[1], q[2]) ≈ gate("cx", q[1], q[2])
    @test gate("cp", q[1], q[2]; cargs=(θ,)) ≈ control(gate("p", q[2]; cargs=(θ,)), q[1])
    @test gate("ccx", q[1], q[2], q[3]) ≈ control(gate("x", q[3]), q[1], q[2])
    @test gate("c3x", q[1], q[2], q[3], q[4]) ≈ control(gate("x", q[4]), q[1], q[2], q[3])
    @test gate("c4x", q[1], q[2], q[3], q[4], q[5]) ≈
        control(gate("x", q[5]), q[1], q[2], q[3], q[4])

    @test gate("cy", q[1], q[2]) ≈ control(gate("y", q[2]), q[1])
    @test gate("cz", q[1], q[2]) ≈ control(gate("z", q[2]), q[1])
    @test gate("ch", q[1], q[2]) ≈ control(gate("h", q[2]), q[1])

    @test gate("crx", q[1], q[2]; cargs=(θ,)) ≈ control(gate("rx", q[2]; cargs=(θ,)), q[1])
    @test gate("cry", q[1], q[2]; cargs=(θ,)) ≈ control(gate("ry", q[2]; cargs=(θ,)), q[1])
    @test gate("crz", q[1], q[2]; cargs=(θ,)) ≈ control(gate("rz", q[2]; cargs=(θ,)), q[1])
    @test gate("crz", q[1], q[2]; cargs=(θ,)) ≈ apply(
        gate("u1", q[2]; cargs=(θ/2,)),
        apply(
            gate("cx", q[1], q[2]),
            apply(gate("u1", q[2]; cargs=(-θ/2,)), gate("cx", q[1], q[2])),
        ),
    )

    @test gate("cswap", q[1], q[2], q[3]) ≈ control(gate("swap", q[2], q[3]), q[1])
    @test gate("cu1", q[1], q[2]; cargs=(θ,)) ≈ control(gate("u1", q[2]; cargs=(θ,)), q[1])
    @test gate("cu3", q[1], q[2]; cargs=(θ, ϕ, λ)) ≈
        control(gate("u3", q[2]; cargs=(θ, ϕ, λ)), q[1])

    return nothing
end

function test_gate_adjoints()
    N = 5
    sq = siteinds("Qubit", N)
    svq = siteinds("vQubit", N)

    for sites in [sq, svq]
        g = gate("cp", sites, 1, 2; cargs=(pi * rand()))
        @test apply(g, gateadjoint(g)) ≈
            ITensors.op("Id", sites, 1) * ITensors.op("Id", sites, 2)

        g = gate("c3x", sites, 2, 3, 4, 5)
        @test apply(g, gateadjoint(g)) ≈
            ITensors.op("Id", sites, 2) *
              ITensors.op("Id", sites, 3) *
              ITensors.op("Id", sites, 4) *
              ITensors.op("Id", sites, 5)

        g = gate("sxdg", sites, 3)
        @test apply(g, gateadjoint(g)) ≈ ITensors.op("Id", sites, 3)

        g = gate("crz", sites, 1, 2; cargs=(pi * rand()))
        @test apply(g, gateadjoint(g)) ≈
            ITensors.op("Id", sites, 1) * ITensors.op("Id", sites, 2)
    end

    return nothing
end

function test_gates_vqubit()
    # Test that the vQubit version of some gates act correctly.
    N = 2
    sites_vec = siteinds("vQubit", N)
    sites = siteinds("Qubit", N)
    v0 = MPS(sites, "0")
    v0_vec = MPS(sites_vec, "0")

    Hv0 = apply(gate("h", sites, 1), v0)
    Hv0_vec = apply(gate("h", sites_vec, 1), v0_vec)  # this is H|0⟩⟨0|H*
    @test Hv0_vec ≈ vec_projector(Hv0; existing_sites=sites_vec)

    CXv0 = apply(gate("cx", sites, 2, 1), v0)
    CXv0_vec = apply(gate("cx", sites_vec, 2, 1), v0_vec)  # this is CX|0⟩⟨0|CX*
    @test CXv0_vec ≈ vec_projector(CXv0; existing_sites=sites_vec)

    return nothing
end

function test_ecr_vqubit()
    # Check that the Qubit and the vQubit version of a gate which is defined from an
    # OpenQASM string match.
    ecr_txt = """OPENQASM 2.0;
  include "qelib1.inc";
  gate rzx(param0) q0,q1 { h q1; cx q0,q1; rz(param0) q1; cx q0,q1; h q1; }
  gate ecr q0,q1 { rzx(pi/7) q0,q1; x q0; rzx(-pi/7) q0,q1; }
  qreg q[4];
  ecr q[0],q[1];
  ecr q[2],q[3];
  ecr q[1],q[2];"""
    s, gs = gates(OpenQASM.parse(ecr_txt), "Qubit"; warn_on_gate_redefinition=false)
    v = MPS(s, "0")
    for g in gs
        v = apply(g, v)
    end
    expvals_qubit = ComplexF64[]
    append!(expvals_qubit, expect(v, "X"))
    append!(expvals_qubit, expect(v, "Y"))
    append!(expvals_qubit, expect(v, "Z"))

    s, gs = gates(OpenQASM.parse(ecr_txt), "vQubit"; warn_on_gate_redefinition=false)
    v = MPS(s, "0")
    for g in gs
        v = apply(g, v)
    end
    expvals_vqubit = ComplexF64[]
    for statename in ["X", "Y", "Z"]
        for i in 1:length(s)
            op_str = repeat(["Id"], length(s))
            op_str[i] = statename
            o = MPS(s, op_str)
            push!(expvals_vqubit, dot(o, v))
        end
    end

    return expvals_qubit ≈ expvals_vqubit
end
