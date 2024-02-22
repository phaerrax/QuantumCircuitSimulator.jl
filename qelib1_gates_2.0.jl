# OpenQASM 2.0 --> ITensor gate map
# ---------------------------------
# 
# The gate set is taken from the Quantum Experience standard header, and adheres to the
# OpenQASM 2.0 specification.

using ITensors

ITensors.op(::OpName"Id", ::SiteType"Qubit") = [
    1 0
    0 1
]
ITensors.op(::OpName"id", st::SiteType"Qubit") = ITensors.op(OpName("Id"), st)

"""
    u_relphase(θ::Real, ϕ::Real, λ::Real)

Different libraries define the three-parameter SU(2) gate differently, i.e. up to a
phase. This function returns the phase ``ℯ^(i f(θ, ϕ, λ))`` by which ITensors' `Rn` operator
and OpenQASM 3.0's `u` gate differ, that is such that
```math
u(θ, ϕ, λ) = u_relphase(θ, ϕ, λ) * Rn(θ, ϕ, λ)
```
"""
function u_relphase(θ::Real, ϕ::Real, λ::Real)
    return exp(-im / 2 * θ)
end

"""
    ITensors.op(::OpName"U", ::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)

Return a generic single-qbit SU(2) gate, as defined by the OpenQASM 3.0 specs:
https://openqasm.com/language/gates.html#built-in-single-qubit-gate-u
"""
function ITensors.op(::OpName"U", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    #
    #               1  ⎛     1+e^(iθ)        -ie^(iλ)(1-e^(iθ))  ⎞
    # U(θ, ϕ, λ) := -  ⎜                                         ⎟
    #               2  ⎝ ie^(iϕ)(1-e^(iθ))  e^(i(ϕ+λ))(1+e^(iθ)) ⎠
    #
    # This gate is already implemented by ITensors, albeit with a different phase.
    # We make this explicit by using the `u_relphase` function above.
    return u_relphase(θ, ϕ, λ) * ITensors.op(OpName("Rn"), st; θ=θ, ϕ=ϕ, λ=λ)
end

function ITensors.op(::OpName"CCCCNOT", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("id"), st)
    not = op(OpName("X"), st)
    return kron(proj0, proj0, proj0, proj0, id) +
           kron(proj0, proj0, proj0, proj1, id) +
           kron(proj0, proj0, proj1, proj0, id) +
           kron(proj0, proj0, proj1, proj1, id) +
           kron(proj0, proj1, proj0, proj0, id) +
           kron(proj0, proj1, proj0, proj1, id) +
           kron(proj0, proj1, proj1, proj0, id) +
           kron(proj0, proj1, proj1, proj1, id) +
           kron(proj1, proj0, proj0, proj0, id) +
           kron(proj1, proj0, proj0, proj1, id) +
           kron(proj1, proj0, proj1, proj0, id) +
           kron(proj1, proj0, proj1, proj1, id) +
           kron(proj1, proj1, proj0, proj0, id) +
           kron(proj1, proj1, proj0, proj1, id) +
           kron(proj1, proj1, proj1, proj0, id) +
           kron(proj1, proj1, proj1, proj1, not)
end

function ITensors.op(::OpName"CH", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("id"), st)
    h = op(OpName("H"), st)
    return kron(proj0, id) + kron(proj1, h)
end

function gate_id(sites::Vector{<:Index}, n::Int)
    return ITensors.op("id", sites, n)
end

function gate_u1(sites::Vector{<:Index}, n::Int, λ::Real)
    return ITensors.op("U", sites, n; θ=0, ϕ=0, λ=λ)
end

function gate_u2(sites::Vector{<:Index}, n::Int, ϕ::Real, λ::Real)
    return ITensors.op("U", sites, n; θ=pi / 2, ϕ=ϕ, λ=λ)
end

function gate_u3(sites::Vector{<:Index}, n::Int, θ::Real, ϕ::Real, λ::Real)
    return ITensors.op("U", sites, n; θ=θ, ϕ=ϕ, λ=λ)
end

function gate_u(sites::Vector{<:Index}, n::Int, θ::Real, ϕ::Real, λ::Real)
    return gate_u3(sites, n, θ, ϕ, λ)
end

function gate_cx(sites::Vector{<:Index}, control::Int, target::Int)
    return ITensors.op("CX", sites, control, target)
end

function gate_x(sites::Vector{<:Index}, n::Int)
    return ITensors.op("X", sites, n)
end

function gate_y(sites::Vector{<:Index}, n::Int)
    return ITensors.op("Y", sites, n)
end

function gate_z(sites::Vector{<:Index}, n::Int)
    return ITensors.op("Z", sites, n)
end

function gate_s(sites::Vector{<:Index}, n::Int)
    return ITensors.op("S", sites, n)
end

function gate_sdg(sites::Vector{<:Index}, n::Int)  # Adjoint of "S"
    # S = Phase(-π/2), where
    #
    #             ⎛ 1       0      ⎞
    # Phase(ϕ) =  ⎜                ⎟
    #             ⎝ 0  exp(im * ϕ) ⎠
    #
    # so S* = Phase(π/2)* = Phase(-π/2)
    return ITensors.op("Phase", sites, n; ϕ=-pi / 2)
end

function gate_h(sites::Vector{<:Index}, n::Int)
    return ITensors.op("H", sites, n)
end

function gate_t(sites::Vector{<:Index}, n::Int)
    return ITensors.op("T", sites, n)
end

function gate_tdg(sites::Vector{<:Index}, n::Int)  # Adjoint of "T"
    # T = Phase(π/4), so T* = Phase(π/4)* = Phase(-π/4)
    return ITensors.op("Phase", sites, n; ϕ=-pi / 4)
end

function gate_ccx(sites::Vector{<:Index}, control1::Int, control2::Int, target::Int)
    return ITensors.op("Toffoli", sites, control1, control2, target)
end

function gate_c3x(
    sites::Vector{<:Index}, control1::Int, control2::Int, control3::Int, target::Int
)  # CCCNOT
    return ITensors.op("CCCNOT", sites, control1, control2, control3, target)
end

function gate_c4x(
    sites::Vector{<:Index},
    control1::Int,
    control2::Int,
    control3::Int,
    control4::Int,
    target::Int,
)  # CCCCNOT
    return ITensors.op("CCCCNOT", sites, control1, control2, control3, control4, target)
end

function gate_rx(sites::Vector{<:Index}, n::Int, θ::Number)
    return ITensors.op("Rx", sites, n; θ=θ)
end

function gate_ry(sites::Vector{<:Index}, n::Int, θ::Number)
    return ITensors.op("Ry", sites, n; θ=θ)
end

function gate_rz(sites::Vector{<:Index}, n::Int, θ::Number)
    return ITensors.op("Rz", sites, n; θ=θ)
end

function gate_cy(sites::Vector{<:Index}, control::Int, target::Int)
    return ITensors.op("CY", sites, control, target)
end

function gate_cz(sites::Vector{<:Index}, control::Int, target::Int)
    return ITensors.op("CZ", sites, control, target)
end

function gate_ch(sites::Vector{<:Index}, control::Int, target::Int)
    return ITensors.op("CH", sites, control, target)
end

function gate_swap(sites::Vector{<:Index}, n1::Int, n2::Int)
    return ITensors.op("Swap", sites, n1, n2)
end

function gate_cswap(sites::Vector{<:Index}, control::Int, target1::Int, target2::Int)
    return ITensors.op("CSwap", sites, control, target1, target2)
end

function gate_crx(sites::Vector{<:Index}, control::Int, target::Int, θ::Number)
    return ITensors.op("CRx", sites, control, target; θ=θ)
end

function gate_cry(sites::Vector{<:Index}, control::Int, target::Int, θ::Number)
    return ITensors.op("CRy", sites, control, target; θ=θ)
end

function gate_crz(sites::Vector{<:Index}, control::Int, target::Int, λ::Number)
    return ITensors.op("CRz", sites, control, target; θ=λ)
    # This is the CRz gate implementation as defined in the qelib1.inc file.
    # It gives an identical result:
    #   return apply(
    #       gate_u1(sites, target, λ / 2),
    #       apply(
    #           gate_cx(sites, control, target),
    #           apply(
    #               gate_u1(sites, target, -λ / 2),
    #               gate_cx(sites, control, target)),
    #       ),
    #   )
end

function ITensors.op(::OpName"CU3", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    u = ITensors.op(OpName("U"), st; θ=θ, ϕ=ϕ, λ=λ)
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("id"), st)
    return kron(proj0, id) + kron(proj1, u)
end

function ITensors.op(::OpName"CU", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    return ITensors.op(OpName("CU3"), st; θ=θ, ϕ=ϕ, λ=λ)
end

function ITensors.op(::OpName"CU1", st::SiteType"Qubit"; λ::Real)
    return ITensors.op(OpName("CU3"), st; θ=0, ϕ=0, λ=λ)
end

function gate_cu1(sites::Vector{<:Index}, control::Int, target::Int, λ::Number)
    return ITensors.op("CU1", sites, control, target; λ=λ)
    # This is the cu1 gate implementation as defined in the qelib1.inc file:
    #
    #   apply(
    #       gate_u1(sites, control, λ / 2),
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

function gate_cu3(
    sites::Vector{<:Index}, control::Int, target::Int, θ::Real, ϕ::Real, λ::Real
)
    return ITensors.op("CU3", sites, control, target; θ=θ, ϕ=ϕ, λ=λ)
    # This is the cu3 gate implementation as defined in the qelib1.inc file:
    # FIXME: it seems like this doesn't return |0⟩⟨0| ⊗ I₂ + |1⟩⟨1| ⊗ U₃(θ,ϕ,λ)...
    #
    #   apply(
    #       gate_u3(sites, target, θ / 2, ϕ, 0),
    #       apply(
    #           gate_cx(sites, control, target),
    #           apply(
    #               gate_u3(sites, target, -θ / 2, 0, -(ϕ + λ) / 2),
    #               apply(
    #                   gate_cx(sites, control, target),
    #                   gate_u1(sites, target, (λ - ϕ) / 2) * gate_u1(sites, control, (λ + ϕ) / 2),
    #               ),
    #           ),
    #       ),
    #   )
end

# I don't recognize these gates...
#
#function gate_c3sqrtx(sites::Vector{<:Index}, n::Int)
#    return ITensors.op("", sites, n)
#end
#function gate_rxx(sites::Vector{<:Index}, n::Int)
#    return ITensors.op("", sites, n)
#end
#
#function gate_rzz(sites::Vector{<:Index}, n::Int)
#    return ITensors.op("", sites, n)
#end
#
#function gate_rccx(sites::Vector{<:Index}, n::Int)
#    return ITensors.op("", sites, n)
#end
#
#function gate_rc3x(sites::Vector{<:Index}, n::Int)
#    return ITensors.op("", sites, n)
#end

function arity(gatename::AbstractString)
    arities = Dict(
        "id" => 1,
        "u1" => 1,
        "u2" => 1,
        "u3" => 1,
        "u" => 1,
        "cx" => 2,
        "x" => 1,
        "y" => 1,
        "z" => 1,
        "s" => 1,
        "sdg" => 1,
        "h" => 1,
        "t" => 1,
        "tdg" => 1,
        "ccx" => 3,
        "c3x" => 4,
        "c4x" => 5,
        "rx" => 1,
        "ry" => 1,
        "rz" => 1,
        "cy" => 2,
        "cz" => 2,
        "ch" => 2,
        "swap" => 2,
        "cswap" => 3,
        "crx" => 2,
        "cry" => 2,
        "crz" => 2,
        "cu1" => 2,
        "cu3" => 2,
        "c3sqrtx" => missing,
        "rxx" => missing,
        "rzz" => missing,
        "rccx" => missing,
        "rc3x" => missing,
    )
    return arities[gatename]
end
