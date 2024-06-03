# Qiskit --> ITensors gate map
# ----------------------------
# 
# The gate set is taken from the Quantum Experience standard header, and adheres to the
# OpenQASM 3.0 specification, except for the U gate, which Qiskit apparently defines
# differently (and we follow them in this file).

# The main interface for the `gate` function is
#   gate(gatename::String, s::Index...; kwargs...)
# defined in src/gate.jl. It tries, in turn,
#   gate(GateName(gatename), st, s...; kwargs...)
#   gate(GateName(gatename), st; kwargs...)
# for every SiteType `st` common to all elements in `s`, and returns the first one that does
# not evaluate to `nothing` (or errors out, if none is found).
#
# This `gate` function is used by the automatic parser `parsegate` in the file
# src/qasm_itensor_parser.jl, which _always_ looks for "classical arguments" (i.e. numbers
# given as parameters) in its instruction, and passes a `cargs` keyword argument to `gate`
# containing them.
# If there are no classical arguments, i.e. if `cargs` is an empty list (usually `Any[]`)
# then a version of `gate` with no keyword arguments is called.
# If there are classical arguments, they keyword argument `cargs` will be used.
# By defining `gate` methods as
#   function gate(::GateName"...", ::SiteType"Qubit", s::Index; cargs)
# the function body can then work with these classical arguments.
# It is expected that `cargs` contain the classical arguments in the very same order as in
# the OpenQASM 3.0 specification of the gate.

function ITensors.op(::OpName"U", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    # This is a new operator for the "Qubit" site type.
    # It returns the generic single-qbit SU(2) gate, as defined in Qiskit's documentation at
    # https://docs.quantum.ibm.com/api/qiskit/qiskit.circuit.library.UGate#ugate
    #
    #                ⎛     cos θ        -e^(iλ) sin θ   ⎞
    # U(θ, ϕ, λ) :=  ⎜                                  ⎟
    #                ⎝ e^(iϕ) sin θ    e^(i(ϕ+λ)) cos θ ⎠
    #
    # This gate is exactly the `Rn` gate implemented by ITensors.
    return ITensors.op(OpName("Rn"), st; θ=θ, ϕ=ϕ, λ=λ)
end

function gate(::GateName"id", ::SiteType"Qubit", s::Index)
    return ITensors.op("Id", s)
end

function gate(::GateName"u1", ::SiteType"Qubit", s::Index; cargs)
    λ::Real = cargs[1]
    return ITensors.op("U", s; θ=0, ϕ=0, λ=λ)
end

function gate(::GateName"u2", ::SiteType"Qubit", s::Index; cargs)
    ϕ::Real = cargs[1]
    λ::Real = cargs[2]
    return ITensors.op("U", s; θ=pi / 2, ϕ=ϕ, λ=λ)
end

function gate(::GateName"u3", ::SiteType"Qubit", s::Index; cargs)
    θ::Real = cargs[1]
    ϕ::Real = cargs[2]
    λ::Real = cargs[3]
    return ITensors.op("U", s; θ=θ, ϕ=ϕ, λ=λ)
end

function gate(::GateName"u", st::SiteType"Qubit", s::Index; cargs)
    return gate(GateName("u3"), st, s; cargs)
end

function gate(::GateName"cx", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CX", control, target)
end

function gate(::GateName"x", ::SiteType"Qubit", s::Index)
    return ITensors.op("X", s)
end

function gate(::GateName"y", ::SiteType"Qubit", s::Index)
    return ITensors.op("Y", s)
end

function gate(::GateName"z", ::SiteType"Qubit", s::Index)
    return ITensors.op("Z", s)
end

function gate(::GateName"p", ::SiteType"Qubit", s::Index; cargs)
    ϕ::Real = cargs[1]
    return ITensors.op("Phase", s; ϕ=ϕ)
end

function gate(::GateName"cp", ::SiteType"Qubit", control::Index, target::Index; cargs)
    ϕ::Real = cargs[1]
    return ITensors.op("CPhase", control, target; ϕ=ϕ)
end

function gate(::GateName"s", ::SiteType"Qubit", s::Index)
    return ITensors.op("S", s)
end

function gate(::GateName"sdg", ::SiteType"Qubit", s::Index)  # Adjoint of "S"
    # S = Phase(-π/2), where
    #
    #             ⎛ 1       0      ⎞
    # Phase(ϕ) =  ⎜                ⎟
    #             ⎝ 0  exp(im * ϕ) ⎠
    #
    # so S* = Phase(π/2)* = Phase(-π/2)
    return ITensors.op("Phase", s; ϕ=-pi / 2)
end

function gate(::GateName"h", ::SiteType"Qubit", s::Index)
    return ITensors.op("H", s)
end

function gate(::GateName"t", ::SiteType"Qubit", s::Index)
    return ITensors.op("T", s)
end

function gate(::GateName"tdg", ::SiteType"Qubit", s::Index)  # Adjoint of "T"
    # T = Phase(π/4), so T* = Phase(π/4)* = Phase(-π/4)
    return ITensors.op("Phase", s; ϕ=-pi / 4)
end

function gate(
    ::GateName"ccx", ::SiteType"Qubit", control1::Index, control2::Index, target::Index
)
    return ITensors.op("Toffoli", control1, control2, target)
end

function gate(
    ::GateName"c3x",
    ::SiteType"Qubit",
    control1::Index,
    control2::Index,
    control3::Index,
    target::Index,
)
    return ITensors.op("CCCNOT", control1, control2, control3, target)
end

function gate(
    ::GateName"c4x",
    ::SiteType"Qubit",
    control1::Index,
    control2::Index,
    control3::Index,
    control4::Index,
    target::Index,
)
    return ITensors.op("CCCCNOT", control1, control2, control3, control4, target)
end

function gate(::GateName"rx", ::SiteType"Qubit", s::Index; cargs)
    θ::Real = cargs[1]
    return ITensors.op("Rx", s; θ=θ)
end

function gate(::GateName"ry", ::SiteType"Qubit", s::Index; cargs)
    θ::Real = cargs[1]
    return ITensors.op("Ry", s; θ=θ)
end

function gate(::GateName"rz", ::SiteType"Qubit", s::Index; cargs)
    θ::Real = cargs[1]
    return ITensors.op("Rz", s; θ=θ)
end

function gate(::GateName"sx", st::SiteType"Qubit", s::Index)
    sdg = ITensors.op(OpName("Phase"), st; ϕ=-pi / 2)
    h = ITensors.op(OpName("H"), st)
    return op(exp(im * pi / 4) * sdg * h * sdg, s)
end

function gate(::GateName"sxdg", st::SiteType"Qubit", s::Index)
    sdg = ITensors.op(OpName("Phase"), st; ϕ=-pi / 2)
    h = ITensors.op(OpName("H"), st)
    return op((exp(im * pi / 4) * sdg * h * sdg)', s)
end

function gate(::GateName"cy", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CY", control, target)
end

function gate(::GateName"cz", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CZ", control, target)
end

function gate(::GateName"ch", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CH", control, target)
end

function gate(::GateName"swap", ::SiteType"Qubit", s1::Index, s2::Index)
    return ITensors.op("Swap", s1, s2)
end

function gate(
    ::GateName"cswap", ::SiteType"Qubit", control::Index, target1::Index, target2::Index
)
    return ITensors.op("CSwap", control, target1, target2)
end

function gate(::GateName"crx", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ::Real = cargs[1]
    return ITensors.op("CRx", control, target; θ=θ)
end

function gate(::GateName"cry", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ::Real = cargs[1]
    return ITensors.op("CRy", control, target; θ=θ)
end

function gate(::GateName"crz", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ::Real = cargs[1]
    return ITensors.op("CRz", control, target; θ=θ)
    # This is the CRz gate implementation as defined in the qelib1.inc file.
    # It gives an identical result:
    #   return apply(
    #       gate_u1(sites, target, λ / 2),
    #       apply(
    #           gate_cx(control, target),
    #           apply(
    #               gate_u1(sites, target, -λ / 2),
    #               gate_cx(control, target)),
    #       ),
    #   )
end

function gate(::GateName"cu1", ::SiteType"Qubit", control::Index, target::Index; cargs)
    λ::Real = cargs[1]
    return ITensors.op("CU1", control, target; λ=λ)
    # This is the cu1 gate implementation as defined in the qelib1.inc file:
    #
    #   apply(
    #       gate_u1(control, λ / 2),
    #       apply(
    #           gate_cx(sites, target, control),
    #           apply(
    #               gate_u1(sites, target, -λ / 2),
    #               apply(
    #                   gate_cx(sites, target, control),
    #                   gate_u1(sites, target, λ / 2)),
    #           ),
    #       ),
    #   )
    #
    # It gives |0⟩⟨0| ⊗ I₂ + |1⟩⟨1| ⊗ U₁(λ) when
    #
    #           ⎛ 1     0    ⎞
    #   U₁(λ) = ⎜            ⎟
    #           ⎝ 0  ℯ^(i λ) ⎠
    #
    # as in the OpenQASM 3.0 specs.
end

function gate(::GateName"cu3", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ::Real = cargs[1]
    ϕ::Real = cargs[2]
    λ::Real = cargs[3]
    return ITensors.op("CU3", control, target; θ=θ, ϕ=ϕ, λ=λ)
    # This is the cu3 gate implementation as defined in the qelib1.inc file:
    # FIXME: it seems like this doesn't return |0⟩⟨0| ⊗ I₂ + |1⟩⟨1| ⊗ U₃(θ,ϕ,λ)...
    #
    #   apply(
    #       gate_u3(target, θ / 2, ϕ, 0),
    #       apply(
    #           gate_cx(control, target),
    #           apply(
    #               gate_u3(target, -θ / 2, 0, -(ϕ + λ) / 2),
    #               apply(
    #                   gate_cx(control, target),
    #                   gate_u1(target, (λ - ϕ) / 2) * gate_u1(control, (λ + ϕ) / 2),
    #               ),
    #           ),
    #       ),
    #   )
end

# I don't recognize these gates...
#
#function gate_c3sqrtx(s::Index)
#    return ITensors.op("", s)
#end
#function gate_rxx(s::Index)
#    return ITensors.op("", s)
#end
#
#function gate_rzz(s::Index)
#    return ITensors.op("", s)
#end
#
#function gate_rccx(s::Index)
#    return ITensors.op("", s)
#end
#
#function gate_rc3x(s::Index)
#    return ITensors.op("", s)
#end
