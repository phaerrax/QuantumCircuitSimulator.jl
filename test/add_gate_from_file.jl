function add_gate_from_file()
    str = "OPENQASM 2.0;
    qreg q[2];
    gate test(a, b, c) q0, q1 {
      h q1;
      cx q0, q1;
      rz(a) q1;
      u2(b, c) q0;
      y q1;
    }"
    decl = OpenQASM.parse(str)
    st = SiteType("Qubit")
    qs = TEM.qbitsites(decl, "Qubit")
    eval(Meta.parse(TEM.definition(decl.prog[2], st)))
    return !isnothing(invokelatest(gate, "test", qs, 1, 2; cargs=(pi / 4, 0, sin(2))))
end
