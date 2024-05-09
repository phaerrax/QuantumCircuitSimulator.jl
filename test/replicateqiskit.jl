const np = pyimport("numpy")
const qiskit = pyimport("qiskit")
const qi = pyimport("qiskit.quantum_info")
const aer = pyimport("qiskit_aer")
const qasm2 = pyimport("qiskit.qasm2")

function qiskitcircuit_noentanglement(output_qasmfile)
    qr = qiskit.QuantumRegister(2, "q")
    qc = qiskit.QuantumCircuit(qr)

    qc.x(qr[0])
    qc.y(qr[1])

    state = qi.Statevector.from_instruction(qc)
    qasm2.dump(qc, output_qasmfile)

    return np.real(state.data), np.imag(state.data)
end

function qiskitcircuit_unitarygates(output_qasmfile)
    sim = aer.UnitarySimulator()
    qr = qiskit.QuantumRegister(1, "q")
    qc = qiskit.QuantumCircuit(qr)

    qc.u(pi / 3, pi / 3, pi / 3, qr[0])

    job = sim.run(qc)
    result = job.result()
    print(result.get_unitary(qc, 4))
    state = qi.Statevector.from_instruction(qc)

    qasm2.dump(qc, output_qasmfile)

    return np.real(state.data), np.imag(state.data)
end

function qiskitcircuit_registermapping(output_qasmfile)
    qr = qiskit.QuantumRegister(5, "q")
    anc = qiskit.QuantumRegister(1, "ancilla")
    qc = qiskit.QuantumCircuit(qr, anc)

    qc.u(pi / 4, 0, pi / 8, qr[1])
    qc.u(pi / 8, pi / 4, 0, qr[2])
    qc.u(0, pi / 8, pi / 4, qr[4])

    qc.u(pi / 2, pi / 2, pi / 6, anc[0])

    qc.x(anc[0])
    qc.h(qr[2])
    qc.t(qr[1])
    qc.cx(qr[1], qr[3])
    qc.cz(qr[0], anc[0])
    qc.cz(qr[1], anc[0])
    qc.cz(qr[2], anc[0])

    qc.ccx(qr[1], qr[2], anc[0])
    qc.h(qr[2])
    qc.h(qr[3])
    qc.swap(qr[1], qr[3])

    state = qi.Statevector.from_instruction(qc)

    qasm2.dump(qc, output_qasmfile)

    return np.real(state.data), np.imag(state.data)
end

function replicateqiskit(qiskitcircuit::Function)
    qasmfile_path, qasmfile_io = mktemp()

    qiskitstate_re, qiskitstate_im = qiskitcircuit(qasmfile_path)
    qiskitstate = complex.(qiskitstate_re, qiskitstate_im)

    # Now we read the OpenQASM file and build the circuit with our library.
    sites, gatelist = open(qasmfile_path, "r") do f
        code = read(f, String)
        gates(code, "Qubit")
    end

    psi = MPS(sites, "0")
    for g in gatelist
        psi = apply(g, psi)
    end
    cbvec = TEM.qiskitvector(psi)

    return isapprox(cbvec, qiskitstate)
end
