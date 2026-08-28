using QuantumCircuitSimulator: instructionsites, freesites

function test_parse_big_circuit()
    circuit_str = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[19];
x q[2];
x q[4];
x q[6];
x q[8];
x q[10];
x q[12];
x q[14];
x q[16];
x q[18];
rz(pi/2) q[0];
rz(pi/2) q[1];
rz(pi/2) q[3];
rz(pi/2) q[5];
rz(pi/2) q[7];
rz(pi/2) q[9];
rz(pi/2) q[11];
rz(pi/2) q[13];
rz(pi/2) q[15];
rz(pi/2) q[17];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
sx q[0];
sx q[1];
sx q[3];
sx q[5];
sx q[7];
sx q[9];
sx q[11];
sx q[13];
sx q[15];
sx q[17];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
rz(pi/2) q[0];
rz(pi/2) q[1];
rz(pi/2) q[3];
rz(pi/2) q[5];
rz(pi/2) q[7];
rz(pi/2) q[9];
rz(pi/2) q[11];
rz(pi/2) q[13];
rz(pi/2) q[15];
rz(pi/2) q[17];
sx q[2];
sx q[4];
sx q[6];
sx q[8];
sx q[10];
sx q[12];
sx q[14];
sx q[16];
sx q[18];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
rz(-pi/2) q[5];
rz(-pi/2) q[7];
rz(-pi/2) q[9];
rz(-pi/2) q[11];
rz(-pi/2) q[13];
rz(-pi/2) q[15];
rz(-pi/2) q[17];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
rz(pi/2) q[1];
rz(pi/2) q[5];
rz(pi/2) q[9];
rz(pi/2) q[11];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
sx q[1];
sx q[5];
sx q[9];
sx q[11];
rz(-pi/2) q[2];
rz(-pi/2) q[6];
rz(-pi/2) q[10];
rz(-pi/2) q[12];
rz(-pi/2) q[1];
rz(-pi/2) q[5];
rz(-pi/2) q[9];
rz(-pi/2) q[11];
sx q[2];
sx q[6];
sx q[10];
sx q[12];
rz(pi/2) q[2];
rz(pi/2) q[6];
rz(pi/2) q[10];
rz(pi/2) q[12];
cx q[3],q[4];
cx q[7],q[8];
cx q[13],q[14];
cx q[15],q[16];
cx q[17],q[18];
cx q[2],q[1];
cx q[6],q[5];
cx q[10],q[9];
cx q[12],q[11];
rz(pi/2) q[11];
rz(pi/2) q[12];
rz(pi/2) q[9];
rz(pi/2) q[10];
rz(pi/2) q[5];
rz(pi/2) q[6];
rz(pi/2) q[1];
rz(pi/2) q[2];
x q[17];
x q[15];
x q[13];
x q[7];
x q[3];
sx q[11];
sx q[12];
sx q[9];
sx q[10];
sx q[5];
sx q[6];
sx q[1];
sx q[2];
rz(pi/2) q[11];
rz(pi/2) q[12];
rz(pi/2) q[9];
rz(pi/2) q[10];
rz(pi/2) q[5];
rz(pi/2) q[6];
rz(pi/2) q[1];
rz(pi/2) q[2];
x q[11];
x q[9];
x q[5];
x q[1];"""
    circ = QuantumCircuit(OpenQASM.parse(circuit_str))
    # Depth has been calculated with `QuantumCircuit` at the time this test function was
    # written. Ideally we should determine it with other means. If anything, this function
    # tests that the circuit gets actually parsed, that a QuantumCircuit object is returned,
    # and that the output (well, the circuit dimensions) have not changed over time as
    # a result of some change in the code.
    @test depth(circ) == 14 && nqbits(circ) == 19
end

function test_parse_instructions()
    circuit_str = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
x q[1];
ccx q[2], q[3], q[4];
cry q[5], q[6];"""
    circ = OpenQASM.parse(circuit_str)
    circ_instructions = filter(line -> line isa OpenQASM.Types.Instruction, circ.prog)
    sites = qbitsites(circ, "Qubit")
    # Note that sites[k] corresponds to qubit q[k-1] in the circuit, not q[k].
    q1 = sites[[2]]
    q234 = sites[[3, 4, 5]]
    q56 = sites[[6, 7]]
    @test instructionsites(circ_instructions[1], sites) == q1
    @test instructionsites(circ_instructions[2], sites) == q234
    @test instructionsites(circ_instructions[3], sites) == q56

    # Robustness check: we throw some fake indices in the mix and check if the result
    # doesn't change (`findsites` should not return the indices based on their position in
    # the list).
    fake_inds = [addtags(siteind("Qubit"), "Site,q[$n]") for n in ['a', 'b']]
    sites = [fake_inds[1]; sites[1:3]; fake_inds[2]; sites[4:end]]
    @test instructionsites(circ_instructions[1], sites) == q1
    @test instructionsites(circ_instructions[2], sites) == q234
    @test instructionsites(circ_instructions[3], sites) == q56
end

function test_barrier_sites()
    circuit_str = """OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[4];
    x q[0];
    barrier q[1], q[3];"""
    circ = OpenQASM.parse(circuit_str)
    barriers = filter(line -> line isa OpenQASM.Types.Barrier, circ.prog)
    sites = qbitsites(circ, "Qubit")
    # Note that sites[k] corresponds to qubit q[k-1] in the circuit, not q[k].
    @test instructionsites(barriers[1], sites) == sites[[2, 4]]
end

function test_parse_barriers()
    # A barrier only constrains the qbits it acts on: gates on unrelated qbits should keep
    # sliding back as usual, right through the barrier.
    circuit_str = """OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[4];
    x q[0];
    x q[2];
    barrier q[0], q[2];
    x q[0];
    x q[1];
    x q[2];"""
    circ = QuantumCircuit(OpenQASM.parse(circuit_str))
    q0, q1, q2, q3 = siteinds(circ)
    @test depth(circ) == 2
    # Layer 1: the first `x q[0]` and `x q[2]`, plus the `x q[1]` that slid all the way
    # back, since it wasn't affected by the barrier.
    @test freesites(circ, circ[1]) == [q3]
    # Layer 2: only the second `x q[0]` and `x q[2]`, stopped by the barrier.
    @test freesites(circ, circ[2]) == [q1, q3]

    # A barrier should slide back as far as possible, like a gate would.  It must only
    # depend on the qbits it actually acts on. Here, the barrier is on q[2], which hasn't
    # been touched yet, so it doesn't force the following `x q[2]` into a new layer even
    # though q[0] and q[1] already have one gate each.
    circuit_str = """OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[3];
    x q[0];
    x q[1];
    barrier q[2];
    x q[2];"""
    circ = QuantumCircuit(OpenQASM.parse(circuit_str))
    @test depth(circ) == 1

    # Nothing can slide past a barrier over every qbit in the circuit.
    circuit_str = """OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    x q[0];
    barrier q[0], q[1];
    x q[0];
    x q[1];"""
    circ = QuantumCircuit(OpenQASM.parse(circuit_str))
    q0, q1 = siteinds(circ)
    @test depth(circ) == 2
    @test freesites(circ, circ[1]) == [q1]
    @test freesites(circ, circ[2]) == ITensors.Index[]
end

function test_empty_gate()
    circ = """OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[3];
    gate empty(a) q1, q2 {}
    x q[1];
    empty(1) q[0], q[2];"""
    circ = OpenQASM.parse(circ_str)
    sites, gate_list = gates(circ, "Qubit"; warn_on_gate_redefinition=true)
    @test gate_list[2] ≈ ITensors.op(I, sites[1], sites[3])
end

function test_layers_outside_domain()
    # Test the inner QuantumCircuit constructor by supplying an array of sites and some
    # gates which act on other sites: it should throw an error.
    circ = """OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[3];
    gate empty(a) q1, q2 {}
    x q[1];
    cx q[0], q[2];
    crz(pi/3) q[1], q[2];
    sx q[0];"""

    # First we build the quantum circuit normally, then we remove one of the sites.
    qc = QuantumCircuit(circ)
    sites = siteinds(qc)

    sites = sites[[1, 3]]
    @test_throws "gates and sites do not match:" QuantumCircuit(sites, instructions(qc))
end
