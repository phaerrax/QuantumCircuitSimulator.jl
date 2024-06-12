using ITensors.SiteTypes: _sitetypes

function ITensors.state(sn::StateName"X+", st::SiteType"vQubit")
    v = ITensors.state(sn, SiteType("Qubit"))
    return LindbladVectorizedTensors.vec(kron(v, v'), LindbladVectorizedTensors.ptmbasis(1))
end

function ITensors.op(::OpName"Bell2Entangler", st::SiteType"Qubit")
    sigma_x = ITensors.op(OpName("X"), st)
    id = ITensors.op(OpName("Id"), st)
    return (1 / sqrt(2)) .* (kron(sigma_x, id) + kron(id, sigma_x))
end

bellstate(::SiteType, sites) = nothing

function bellstate(sites)
    common_stypes = _sitetypes(ITensors.commontags(sites...))
    for st in common_stypes
        v = bellstate(st, sites)
        !isnothing(v) && return v
    end
    return throw(
        ArgumentError(
            "Overload of \"kicked_ising_gate_seq\" function not found for " *
            "Index tags $(ITensors.tags.(sites))",
        ),
    )
end

function bellstate(::SiteType"Qubit", sites)
    N = length(sites)
    v = MPS(sites, [j == 1 ? "X+" : "0" for j in 1:N])
    # Apply the "Bell entangler" gate pair by pair to create the initial state.
    for j in 2:2:N
        v = apply(op("Bell2Entangler", sites, j, j + 1), v)
    end
    return v
end

function bellstate(::SiteType"vQubit", sites)
    N = length(sites)
    v = MPS(sites, [j == 1 ? "X+" : "0" for j in 1:N])
    # Apply the "Bell entangler" gate pair by pair to create the initial state.
    for j in 2:2:N
        entangler = LindbladVectorizedTensors.adjointmap_itensor(
            "Bell2Entangler", sites, j, j + 1
        )
        v = apply(entangler, v)
    end
    return v
end

function bellstate_openqasm(; atol=eps(Float64))
    init0 = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[37];
x q[2];
x q[4];
x q[6];
x q[8];
x q[10];
x q[12];
x q[14];
x q[16];
x q[18];
x q[20];
x q[22];
x q[24];
x q[26];
x q[28];
x q[30];
x q[32];
x q[34];
x q[36];
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
rz(pi/2) q[19];
rz(pi/2) q[21];
rz(pi/2) q[23];
rz(pi/2) q[25];
rz(pi/2) q[27];
rz(pi/2) q[29];
rz(pi/2) q[31];
rz(pi/2) q[33];
rz(pi/2) q[35];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
rz(pi/2) q[20];
rz(pi/2) q[22];
rz(pi/2) q[24];
rz(pi/2) q[26];
rz(pi/2) q[28];
rz(pi/2) q[30];
rz(pi/2) q[32];
rz(pi/2) q[34];
rz(pi/2) q[36];
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
sx q[19];
sx q[21];
sx q[23];
sx q[25];
sx q[27];
sx q[29];
sx q[31];
sx q[33];
sx q[35];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
rz(pi/2) q[20];
rz(pi/2) q[22];
rz(pi/2) q[24];
rz(pi/2) q[26];
rz(pi/2) q[28];
rz(pi/2) q[30];
rz(pi/2) q[32];
rz(pi/2) q[34];
rz(pi/2) q[36];
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
rz(pi/2) q[19];
rz(pi/2) q[21];
rz(pi/2) q[23];
rz(pi/2) q[25];
rz(pi/2) q[27];
rz(pi/2) q[29];
rz(pi/2) q[31];
rz(pi/2) q[33];
rz(pi/2) q[35];
sx q[2];
sx q[4];
sx q[6];
sx q[8];
sx q[10];
sx q[12];
sx q[14];
sx q[16];
sx q[18];
sx q[20];
sx q[22];
sx q[24];
sx q[26];
sx q[28];
sx q[30];
sx q[32];
sx q[34];
sx q[36];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
rz(-pi/2) q[5];
rz(-pi/2) q[7];
rz(-pi/2) q[9];
rz(-pi/2) q[11];
rz(-pi/2) q[13];
rz(-pi/2) q[15];
rz(-pi/2) q[17];
rz(-pi/2) q[19];
rz(-pi/2) q[21];
rz(-pi/2) q[23];
rz(-pi/2) q[25];
rz(-pi/2) q[27];
rz(-pi/2) q[29];
rz(-pi/2) q[31];
rz(-pi/2) q[33];
rz(-pi/2) q[35];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
rz(pi/2) q[20];
rz(pi/2) q[22];
rz(pi/2) q[24];
rz(pi/2) q[26];
rz(pi/2) q[28];
rz(pi/2) q[30];
rz(pi/2) q[32];
rz(pi/2) q[34];
rz(pi/2) q[36];
rz(pi/2) q[1];
rz(pi/2) q[5];
rz(pi/2) q[9];
rz(pi/2) q[11];
rz(pi/2) q[21];
rz(pi/2) q[33];
rz(pi/2) q[35];
rz(pi/2) q[2];
rz(pi/2) q[4];
rz(pi/2) q[6];
rz(pi/2) q[8];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[14];
rz(pi/2) q[16];
rz(pi/2) q[18];
rz(pi/2) q[20];
rz(pi/2) q[22];
rz(pi/2) q[24];
rz(pi/2) q[26];
rz(pi/2) q[28];
rz(pi/2) q[30];
rz(pi/2) q[32];
rz(pi/2) q[34];
rz(pi/2) q[36];
sx q[1];
sx q[5];
sx q[9];
sx q[11];
sx q[21];
sx q[33];
sx q[35];
rz(-pi/2) q[2];
rz(-pi/2) q[6];
rz(-pi/2) q[10];
rz(-pi/2) q[12];
rz(-pi/2) q[22];
rz(-pi/2) q[34];
rz(-pi/2) q[36];
rz(-pi/2) q[1];
rz(-pi/2) q[5];
rz(-pi/2) q[9];
rz(-pi/2) q[11];
rz(-pi/2) q[21];
rz(-pi/2) q[33];
rz(-pi/2) q[35];
sx q[2];
sx q[6];
sx q[10];
sx q[12];
sx q[22];
sx q[34];
sx q[36];
rz(pi/2) q[2];
rz(pi/2) q[6];
rz(pi/2) q[10];
rz(pi/2) q[12];
rz(pi/2) q[22];
rz(pi/2) q[34];
rz(pi/2) q[36];"""
    init1_ecr = """OPENQASM 2.0;
include "qelib1.inc";
gate rzx(param0) q0,q1 { h q1; cx q0,q1; rz(param0) q1; cx q0,q1; h q1; }
gate ecr q0,q1 { rzx(pi/4) q0,q1; x q0; rzx(-pi/4) q0,q1; }
qreg q[37];
ecr q[3],q[4];
ecr q[7],q[8];
ecr q[13],q[14];
ecr q[15],q[16];
ecr q[17],q[18];
ecr q[19],q[20];
ecr q[23],q[24];
ecr q[25],q[26];
ecr q[27],q[28];
ecr q[29],q[30];
ecr q[31],q[32];
ecr q[2],q[1];
ecr q[6],q[5];
ecr q[10],q[9];
ecr q[12],q[11];
ecr q[22],q[21];
ecr q[34],q[33];
ecr q[36],q[35];"""
    init2 = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[37];
rz(pi/2) q[35];
rz(pi/2) q[36];
rz(pi/2) q[33];
rz(pi/2) q[34];
rz(pi/2) q[21];
rz(pi/2) q[22];
rz(pi/2) q[11];
rz(pi/2) q[12];
rz(pi/2) q[9];
rz(pi/2) q[10];
rz(pi/2) q[5];
rz(pi/2) q[6];
rz(pi/2) q[1];
rz(pi/2) q[2];
x q[31];
x q[29];
x q[27];
x q[25];
x q[23];
x q[19];
x q[17];
x q[15];
x q[13];
x q[7];
x q[3];
sx q[35];
sx q[36];
sx q[33];
sx q[34];
sx q[21];
sx q[22];
sx q[11];
sx q[12];
sx q[9];
sx q[10];
sx q[5];
sx q[6];
sx q[1];
sx q[2];
rz(pi/2) q[35];
rz(pi/2) q[36];
rz(pi/2) q[33];
rz(pi/2) q[34];
rz(pi/2) q[21];
rz(pi/2) q[22];
rz(pi/2) q[11];
rz(pi/2) q[12];
rz(pi/2) q[9];
rz(pi/2) q[10];
rz(pi/2) q[5];
rz(pi/2) q[6];
rz(pi/2) q[1];
rz(pi/2) q[2];
x q[35];
x q[33];
x q[21];
x q[11];
x q[9];
x q[5];
x q[1];"""
    sites, qmap, g0 = gates(OpenQASM.parse(init0), "Qubit")
    nsites = length(sites)
    g1 = gates(OpenQASM.parse(init1_ecr), sites, qmap)
    g2 = gates(OpenQASM.parse(init2), sites, qmap)
    openqasm_plusbell = MPS(sites, "0")
    for g in [g0; g1; g2]
        openqasm_plusbell = apply(g, openqasm_plusbell)
    end

    itensor_plusbell = bellstate(sites)

    p_itensor = projector(itensor_plusbell)
    p_openqasm = projector(openqasm_plusbell)
    return norm(p_itensor - p_openqasm) < atol
end
