using TEM, ITensors

let
    # Qiskit's way of writing the result.
    str_to_c(s) = complex(parse.(Float64, split(s, ","))...)
    v = open("test/qiskit_basis_ordering_final_state.dat") do f
        return str_to_c.(eachline(f))
    end

    # Now we read the OpenQASM file and build the circuit with our library.
    circuitfile = "test/qiskit_basis_ordering_circuit.qasm2"
    sites, gatelist = open(circuitfile, "r") do f
        code = read(f, String)
        gates(code, "Qubit")
    end

    psi = MPS(sites, "0")
    for g in gatelist
        psi = apply(g, psi)
    end
    cbvec = TEM.qiskitvector(psi)
    return isapprox(cbvec, v)
end
