# OpenQASM 3.0 --> ITensor gate map
# ---------------------------------
# 
# The gate set is taken from the Quantum Experience standard header, and adheres to the
# OpenQASM 3.0 specification.

# Define new gates as:
#   gate(::GateName, ::SiteType, ::Index...; kwargs...) -> ITensor
#   gate(::GateName, ::SiteType; kwargs...) -> Julia matrix

"""
    u_relphase_openqasm3(θ::Real, ϕ::Real, λ::Real)

Different libraries define the three-parameter SU(2) gate differently, i.e. up to a
phase. This function returns the phase ``ℯ^(i f(θ, ϕ, λ))`` by which ITensors' `Rn` operator
and OpenQASM 2.0's `u` gate differ, that is such that
```math
u(θ, ϕ, λ) = u_relphase_openqasm3(θ, ϕ, λ) * Rn(θ, ϕ, λ)
```
"""
function u_relphase_openqasm3(θ::Real, ϕ::Real, λ::Real)
    return exp(im / 2 * θ)
end

function ITensors.op(::OpName"U", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    # This is a new operator for the "Qubit" site type.
    # It returns the generic single-qbit SU(2) gate, as defined in the OpenQASM 3.0
    # specification (https://openqasm.com/language/gates.html#U).
    #
    #               1  ⎛     1+e^(iθ)        -ie^(iλ)(1-e^(iθ))  ⎞
    # U(θ, ϕ, λ) := -  ⎜                                         ⎟
    #               2  ⎝ ie^(iϕ)(1-e^(iθ))  e^(i(ϕ+λ))(1+e^(iθ)) ⎠
    #
    # This gate is implemented by ITensors, albeit with a different phase.
    # We make this explicit by using the `u_relphase...` function above.
    return u_relphase_openqasm3(θ, ϕ, λ) * ITensors.op(OpName("Rn"), st; θ=θ, ϕ=ϕ, λ=λ)
end

function gate(::GateName"id", ::SiteType"Qubit", s::Index)
    return ITensors.op("Id", s)
end

function gate(::GateName"u1", ::SiteType"Qubit", s::Index, λ::Real)
    return ITensors.op("U", s; θ=0, ϕ=0, λ=λ)
end

function gate(::GateName"u2", ::SiteType"Qubit", s::Index, ϕ::Real, λ::Real)
    return ITensors.op("U", s; θ=pi / 2, ϕ=ϕ, λ=λ)
end

function gate(::GateName"u3", ::SiteType"Qubit", s::Index, θ::Real, ϕ::Real, λ::Real)
    return ITensors.op("U", s; θ=θ, ϕ=ϕ, λ=λ)
end

function gate(::GateName"u", ::SiteType"Qubit", s::Index, θ::Real, ϕ::Real, λ::Real)
    return gate_u3(s, θ, ϕ, λ)
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

function gate(::GateName"p", ::SiteType"Qubit", s::Index, ϕ::Real)
    return ITensors.op("Phase", s; ϕ=ϕ)
end

function gate(::GateName"cp", ::SiteType"Qubit", control::Index, target::Index, ϕ::Real)
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

function gate(::GateName"rx", ::SiteType"Qubit", s::Index, θ::Real)
    return ITensors.op("Rx", s; θ=θ)
end

function gate(::GateName"ry", ::SiteType"Qubit", s::Index, θ::Real)
    return ITensors.op("Ry", s; θ=θ)
end

function gate(::GateName"rz", ::SiteType"Qubit", s::Index, θ::Real)
    return ITensors.op("Rz", s; θ=θ)
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

function gate(::GateName"crx", ::SiteType"Qubit", control::Index, target::Index, θ::Real)
    return ITensors.op("CRx", control, target; θ=θ)
end

function gate(::GateName"cry", ::SiteType"Qubit", control::Index, target::Index, θ::Real)
    return ITensors.op("CRy", control, target; θ=θ)
end

function gate(::GateName"crz", ::SiteType"Qubit", control::Index, target::Index, λ::Real)
    return ITensors.op("CRz", control, target; θ=λ)
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

function gate(::GateName"cu1", ::SiteType"Qubit", control::Index, target::Index, λ::Real)
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

function gate(
    ::GateName"cu3",
    ::SiteType"Qubit",
    control::Index,
    target::Index,
    θ::Real,
    ϕ::Real,
    λ::Real,
)
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
