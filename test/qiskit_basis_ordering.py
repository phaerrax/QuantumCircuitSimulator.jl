from qiskit import QuantumCircuit, QuantumRegister, qasm2, qasm3
import qiskit.quantum_info as qi
from math import pi

qr = QuantumRegister(3, 'q')
qc = QuantumCircuit(qr)

qc.x(qr[0])
qc.y(qr[1])
qc.h(qr[2])

print(qc.draw())

state = qi.Statevector.from_instruction(qc)

with open("qiskit_basis_ordering_final_state.dat", "w") as out:
    for x in state.data:
        print(f"{x.real},{x.imag}", file=out)

with open("qiskit_basis_ordering_circuit.qasm2", "w") as out:
    qasm2.dump(qc, out)

with open("qiskit_basis_ordering_circuit.qasm3", "w") as out:
    qasm3.dump(qc, out)
