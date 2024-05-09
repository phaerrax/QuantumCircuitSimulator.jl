from qiskit import QuantumCircuit, QuantumRegister, qasm2
import qiskit.quantum_info as qi
from math import pi

qr = QuantumRegister(5, "q")
anc = QuantumRegister(1, "ancilla")
qc = QuantumCircuit(qr, anc)

qc.u(pi / 4, 0, pi / 8, qr[1])
qc.u(pi / 8, pi / 4, 0, qr[2])
qc.u(0, pi / 8, pi / 4, qr[4])

qc.u(pi / 2, pi / 2, pi / 6, anc[0])

qc.x(anc[0])
qc.h(qr[2])
qc.t(qr[1])
qc.cx(qr[1], qr[3])
qc.cz(qr[0:3], anc[0])

qc.ccx(qr[1], qr[2], anc[0])
qc.h(qr[2:4])
qc.swap(qr[1], qr[3])

print(qc.draw())

state = qi.Statevector.from_instruction(qc)

with open("test/sample_circuit_final_state.dat", "w") as out:
    for x in state.data:
        print(f"{x.real},{x.imag}", file=out)

with open("test/sample_circuit.qasm2", "w") as out:
    qasm2.dump(qc, out)
