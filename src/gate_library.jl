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

# Different libraries define the three-parameter SU(2) gate differently, i.e. up to a phase.
# This function returns the phase ``ℯ^(i f(θ, ϕ, λ))`` by which ITensors' `Rn` operator and
# the `u` gate from OpenQASM 2.0/OpenQASM 3.0/Qiskit `u` differ, that is, such that
#   u(θ, ϕ, λ) = phase(θ, ϕ, λ) * Rn(θ, ϕ, λ)
openqasm2_phase(θ, ϕ, λ) = exp(-im / 2 * (ϕ + λ))
openqasm3_phase(θ, ϕ, λ) = exp(im / 2 * θ)
qiskit_phase(θ, ϕ, λ) = 1

# Interface for vectorised states
# -------------------------------

# This retrieves a gate G for the Qubit site type, and then returns the vectorisation of
# the G ⋅ G† operation, that can be applied to (vectorised) mixed states or operators (via
# its adjoint).
function gate(gn::GateName, vst::SiteType"vQubit", vs::Index...; cargs...)
    s = [siteind("Qubit") for _ in vs]
    g = gate(gn, SiteType("Qubit"), s...; cargs...)
    return adjointmap_itensor(g; orig_sites=s, vec_sites=[vs...])
    # Note that both `orig_sites` and `vec_sites` must be _vectors_ of indices: `s` already
    # is so, but we must wrap `vs...` in a vector for this reason.
end

# Basic unitary gates
# -------------------

function gate(::GateName"id", ::SiteType"Qubit", s::Index...)
    return ITensors.op("Id", s...)
end

function ITensors.op(::OpName"U", ::SiteType"Qubit", s::Index; θ::Real, ϕ::Real, λ::Real)
    # This is a new operator for the "Qubit" site type.
    # It returns the generic single-qbit SU(2) gate, as defined in Qiskit's documentation at
    # https://docs.quantum.ibm.com/api/qiskit/qiskit.circuit.library.UGate#ugate
    #
    #                ⎛     cos θ/2        -e^(iλ) sin θ/2  ⎞
    # U(θ, ϕ, λ) :=  ⎜                                     ⎟
    #                ⎝ e^(iϕ) sin θ/2   e^(i(ϕ+λ)) cos θ/2 ⎠
    #
    # This coincides with ITensors' Rn operator.
    return qiskit_phase(θ, ϕ, λ) * ITensors.op("Rn", s; θ, ϕ, λ)
end

# We define the "u" gate through the ITensor "U" operator, then all other parametric unitary
# gates in terms of this gate.
function gate(::GateName"u", ::SiteType"Qubit", s::Index; cargs)
    θ, ϕ, λ=cargs
    return ITensors.op("U", s; θ, ϕ, λ)
end

function gate(::GateName"u1", ::SiteType"Qubit", s::Index; cargs)
    λ = only(cargs)
    return gate("u", s; cargs=(0, 0, λ))
end

function gate(::GateName"u2", ::SiteType"Qubit", s::Index; cargs)
    ϕ, λ = cargs
    return gate("u", s; cargs=(pi/2, ϕ, λ))
end

function gate(::GateName"u3", ::SiteType"Qubit", s::Index; cargs)
    return gate("u", s; cargs)
end

# Basic gates
# -----------

function gate(::GateName"x", ::SiteType"Qubit", s::Index)
    return ITensors.op("X", s)
end

function gate(::GateName"y", ::SiteType"Qubit", s::Index)
    return ITensors.op("Y", s)
end

function gate(::GateName"z", ::SiteType"Qubit", s::Index)
    return ITensors.op("Z", s)
end

function gate(::GateName"h", ::SiteType"Qubit", s::Index)
    return ITensors.op("H", s)
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
    return ITensors.op("Phase", s; ϕ=(-pi / 2))
end

function gate(::GateName"p", ::SiteType"Qubit", s::Index; cargs)
    return ITensors.op("Phase", s; ϕ=only(cargs))
end

function gate(::GateName"t", ::SiteType"Qubit", s::Index)
    return ITensors.op("T", s)
end

function gate(::GateName"tdg", ::SiteType"Qubit", s::Index)  # Adjoint of "T"
    # T = Phase(π/4), so T* = Phase(π/4)* = Phase(-π/4)
    return ITensors.op("Phase", s; ϕ=(-pi / 4))
end

# Rotation gates
# --------------

function gate(::GateName"rx", ::SiteType"Qubit", s::Index; cargs)
    θ = only(cargs)
    return ITensors.op("Rx", s; θ)
end

function gate(::GateName"ry", ::SiteType"Qubit", s::Index; cargs)
    θ = only(cargs)
    return ITensors.op("Ry", s; θ)
end

function gate(::GateName"rz", ::SiteType"Qubit", s::Index; cargs)
    θ = only(cargs)
    return ITensors.op("Rz", s; θ)
end

function gate(::GateName"rxx", ::SiteType"Qubit", q1::Index, q2::Index; cargs)
    # The definition should match
    #      ┌───┐                   ┌───┐
    # q_0: ┤ H ├──■─────────────■──┤ H ├
    #      ├───┤┌─┴─┐┌───────┐┌─┴─┐├───┤
    # q_1: ┤ H ├┤ X ├┤ Rz(0) ├┤ X ├┤ H ├
    #      └───┘└───┘└───────┘└───┘└───┘
    θ = only(cargs)
    rxx = [
        cos(θ/2) 0 0 -im * sin(θ/2)
        0 cos(θ/2) -im * sin(θ/2) 0
        0 -im * sin(θ/2) cos(θ/2) 0
        -im * sin(θ/2) 0 0 cos(θ/2)
    ]
    return gate(rxx, q2, q1)
end

function gate(::GateName"ryy", ::SiteType"Qubit", q1::Index, q2::Index; cargs)
    # The definition should match
    #      ┌─────┐                   ┌────┐
    # q_0: ┤ √X† ├──■─────────────■──┤ √X ├
    #      ├─────┤┌─┴─┐┌───────┐┌─┴─┐├────┤
    # q_1: ┤ √X† ├┤ X ├┤ Rz(θ) ├┤ X ├┤ √X ├
    #      └─────┘└───┘└───────┘└───┘└────┘
    θ = only(cargs)
    ryy = [
        cos(θ/2) 0 0 im * sin(θ/2)
        0 cos(θ/2) -im * sin(θ/2) 0
        0 -im * sin(θ/2) cos(θ/2) 0
        im * sin(θ/2) 0 0 cos(θ/2)
    ]
    return gate(ryy, q2, q1)
end

function gate(::GateName"rzz", ::SiteType"Qubit", q1::Index, q2::Index; cargs)
    # The definition should match
    # q_0: ──■─────────────■──
    #      ┌─┴─┐┌───────┐┌─┴─┐
    # q_1: ┤ X ├┤ Rz(0) ├┤ X ├
    #      └───┘└───────┘└───┘
    θ = only(cargs)
    rzz = diagm(0 => [cis(-θ/2), cis(θ/2), cis(θ/2), cis(-θ/2)])
    return gate(rzz, q2, q1)
end

# Controlled gates
# ----------------

function gate(::GateName"cx", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CX", control, target)
end

function gate(::GateName"cnot", st::SiteType"Qubit", control::Index, target::Index)
    # "cx" alias
    return gate(GateName("cx"), st, control, target)
end

function gate(::GateName"cp", ::SiteType"Qubit", control::Index, target::Index; cargs)
    ϕ = only(cargs)
    return ITensors.op("CPhase", control, target; ϕ=ϕ)
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

function gate(::GateName"cy", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CY", control, target)
end

function gate(::GateName"cz", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CZ", control, target)
end

function gate(::GateName"ch", ::SiteType"Qubit", control::Index, target::Index)
    return ITensors.op("CH", control, target)
end

function gate(::GateName"crx", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ = only(cargs)
    return ITensors.op("CRx", control, target; θ=θ)
end

function gate(::GateName"cry", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ = only(cargs)
    return ITensors.op("CRy", control, target; θ=θ)
end

function gate(::GateName"crz", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ = only(cargs)
    return ITensors.op("CRz", control, target; θ=θ)
    # This is the CRz gate implementation as defined in the qelib1.inc file.
    # It gives an identical result:
    #   return apply(
    #       gate("u1", target; cargs=(λ/2)),
    #       apply(
    #           gate("cx", control, target),
    #           apply(
    #               gate("u1", target; cargs=(-λ/2)),
    #               gate("cx", control, target)),
    #       ),
    #   )
end

function gate(
    ::GateName"cswap", ::SiteType"Qubit", control::Index, target1::Index, target2::Index
)
    return ITensors.op("CSwap", control, target1, target2)
end

function gate(::GateName"cu1", ::SiteType"Qubit", control::Index, target::Index; cargs)
    λ = only(cargs)
    # CU₁(λ) = |0⟩⟨0| ⊗ I₂ + |1⟩⟨1| ⊗ diag(1, e^(iλ)) = diag(1, 1, 1, e^(iλ))
    # which coincides with CPhase.
    return ITensors.op("CPhase", control, target; ϕ=λ)
end

function gate(::GateName"cu3", ::SiteType"Qubit", control::Index, target::Index; cargs)
    θ, ϕ, λ = cargs
    # We build |0⟩⟨0| ⊗ I₂ + |1⟩⟨1| ⊗ U₃(θ,ϕ,λ) directly.
    u = [
        cos(θ / 2) -cis(λ) * sin(θ / 2)
        cis(ϕ) * sin(θ / 2) cis(ϕ + λ) * cos(θ / 2)
    ]
    cu3 = zeros(ComplexF64, 4, 4)
    cu3[1, 1] = 1
    cu3[2, 2] = 1
    cu3[3:4, 3:4] = u
    return gate(cu3, target, control)
end

# Other gates
# -----------

function gate(::GateName"sx", ::SiteType"Qubit", s::Index)
    return ITensors.op("√X", s)
end

function ITensors.op(::OpName"√X†", st::SiteType"Qubit")
    return ITensors.op(OpName("√X"), st)'
end

function gate(::GateName"sxdg", ::SiteType"Qubit", s::Index)
    return ITensors.op("√X†", s)
end

function gate(::GateName"swap", ::SiteType"Qubit", s1::Index, s2::Index)
    return ITensors.op("Swap", s1, s2)
end

function gate(
    ::GateName"rccx", ::SiteType"Qubit", control1::Index, control2::Index, target::Index
)
    # Simplified Toffoli gate.
    # The definition should match
    # q_1: ───────────────────────■───────────────────────
    #                             │
    # q_2: ────────────■──────────┼─────────■─────────────
    #      ┌───┐┌───┐┌─┴─┐┌────┐┌─┴─┐┌───┐┌─┴─┐┌────┐┌───┐
    # q_3: ┤ H ├┤ T ├┤ X ├┤ T† ├┤ X ├┤ T ├┤ X ├┤ T† ├┤ H ├
    #      └───┘└───┘└───┘└────┘└───┘└───┘└───┘└────┘└───┘
    rccx = [
        1 0 0 0 0 0 0 0
        0 1 0 0 0 0 0 0
        0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 -im
        0 0 0 0 1 0 0 0
        0 0 0 0 0 -1 0 0
        0 0 0 0 0 0 1 0
        0 0 0 im 0 0 0 0
    ]
    return gate(rccx, control1, control2, target)
end

function gate(
    ::GateName"rc3x",
    ::SiteType"Qubit",
    control1::Index,
    control2::Index,
    control3::Index,
    target::Index,
)
    # Simplified 3-controlled Toffoli gate.
    # The definition should match
    # q_1: ──────────────────■────────────■────────────────────────────
    #                        │            │
    # q_2: ──────────────────┼─────■──────┼─────■──────────────────────
    #                        │     │      │     │
    # q_3: ────────■─────────┼─────┼──────┼─────┼────────────■─────────
    #       ┌─┐┌─┐┌┴┐┌──┐┌─┐┌┴┐┌─┐┌┴┐┌──┐┌┴┐┌─┐┌┴┐┌──┐┌─┐┌─┐┌┴┐┌──┐┌─┐
    # q_4: ─┤H├┤T├┤X├┤T†├┤H├┤X├┤T├┤X├┤T†├┤X├┤T├┤X├┤T†├┤H├┤T├┤X├┤T†├┤H├─
    #       └─┘└─┘└─┘└──┘└─┘└─┘└─┘└─┘└──┘└─┘└─┘└─┘└──┘└─┘└─┘└─┘└──┘└─┘
    # (https://arxiv.org/abs/1508.03273, Figure 4)
    rc3x = [
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 im 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 -im 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0
    ]
    return gate(rc3x, control1, control2, control3, target)
end
