using LindbladVectorizedTensors: vop

gateadjoint(g) = swapprime(dag(g), 0 => 1)

function gateadjoints()
    tests = Bool[]

    N = 5
    sq = siteinds("Qubit", N)
    svq = siteinds("vQubit", N)

    for sites in [sq, svq]
        g = gate("cp", sites, 1, 2; cargs=(pi * rand()))
        push!(
            tests,
            isapprox(apply(g, gateadjoint(g)), op("Id", sites, 1) * op("Id", sites, 2)),
        )
        g = gate("c3x", sites, 2, 3, 4, 5)
        push!(
            tests,
            isapprox(
                apply(g, gateadjoint(g)),
                op("Id", sites, 2) *
                op("Id", sites, 3) *
                op("Id", sites, 4) *
                op("Id", sites, 5),
            ),
        )
        g = gate("sxdg", sites, 3)
        push!(tests, isapprox(apply(g, gateadjoint(g)), op("Id", sites, 3)))
        g = gate("crz", sites, 1, 2; cargs=(pi * rand()))
        push!(
            tests,
            isapprox(apply(g, gateadjoint(g)), op("Id", sites, 1) * op("Id", sites, 2)),
        )
    end

    return all(tests)
end

function gates_vqubit()
    tests = Bool[]
    # Here we test that the vQubit version of some gates act correctly.
    # H|0⟩ = 1/√2 (|0⟩ + |1⟩)
    N = 2
    sites_vec = siteinds("vQubit", N)
    sites = siteinds("Qubit", N)
    v0 = MPS(sites, "0")
    v0_vec = MPS(sites_vec, "0")

    Hv0 = apply(gate("h", sites, 1), v0)
    Hv0_vec = apply(gate("h", sites_vec, 1), v0_vec)  # this is H|0⟩⟨0|H*
    push!(
        tests,
        isapprox(
            Hv0_vec,
            vec_purestate_densitymatrix(SiteType("Qubit"), Hv0; existing_sites=sites_vec),
        ),
    )

    CXv0 = apply(gate("cx", sites, 2, 1), v0)
    CXv0_vec = apply(gate("cx", sites_vec, 2, 1), v0_vec)  # this is CX|0⟩⟨0|CX*
    push!(
        tests,
        isapprox(
            CXv0_vec,
            vec_purestate_densitymatrix(SiteType("Qubit"), CXv0; existing_sites=sites_vec),
        ),
    )

    return all(tests)
end

function ecr_vqubit()
    # Here we check that the Qubit and the vQubit version of a gate which is defined from
    # an OpenQASM string match.
    ecr_txt = """OPENQASM 2.0;
  include "qelib1.inc";
  gate rzx(param0) q0,q1 { h q1; cx q0,q1; rz(param0) q1; cx q0,q1; h q1; }
  gate ecr q0,q1 { rzx(pi/7) q0,q1; x q0; rzx(-pi/7) q0,q1; }
  qreg q[4];
  ecr q[0],q[1];
  ecr q[2],q[3];
  ecr q[1],q[2];"""
    s, _, gs = gates(OpenQASM.parse(ecr_txt), "Qubit")
    v = MPS(s, "0")
    for g in gs
        v = apply(g, v)
    end
    expvals_qubit = ComplexF64[]
    append!(expvals_qubit, expect(v, "X"))
    append!(expvals_qubit, expect(v, "Y"))
    append!(expvals_qubit, expect(v, "Z"))

    s, _, gs = gates(OpenQASM.parse(ecr_txt), "vQubit")
    v = MPS(s, "0")
    for g in gs
        v = apply(g, v)
    end
    expvals_vqubit = ComplexF64[]
    for statename in ["vX", "vY", "vZ"]
        for i in 1:length(s)
            op_str = repeat(["vId"], length(s))
            op_str[i] = statename
            o = MPS(s, op_str)
            push!(expvals_vqubit, dot(o, v))
        end
    end

    return isapprox(expvals_qubit, expvals_vqubit)
end
