# OpenQASM 3.0 --> ITensor gate map (with vectorized operators)
# -------------------------------------------------------------
# 
# The gate set is taken from the Quantum Experience standard header, and adheres to the
# OpenQASM 3.0 specification.

using ITensors

include("qelib1_gates_3.0.jl")
include("vQubit.jl")

function gate_id(sites::Vector{<:Index}, n::Int)
    return ITensors.op("Id", sites, n)
end

function gate_u1(sites::Vector{<:Index}, n::Int, λ::Real)
    return adjointmap_itensor("U", sites, n; θ=0, ϕ=0, λ=λ)
end

function gate_u2(sites::Vector{<:Index}, n::Int, ϕ::Real, λ::Real)
    return adjointmap_itensor("U", sites, n; θ=pi / 2, ϕ=ϕ, λ=λ)
end

function gate_u3(sites::Vector{<:Index}, n::Int, θ::Real, ϕ::Real, λ::Real)
    return adjointmap_itensor("U", sites, n; θ=θ, ϕ=ϕ, λ=λ)
end

function gate_u(sites::Vector{<:Index}, n::Int, θ::Real, ϕ::Real, λ::Real)
    return gate_u3(sites, n, θ, ϕ, λ)
end

function gate_cx(sites::Vector{<:Index}, control::Int, target::Int)
    return adjointmap_itensor("CX", sites, control, target)
end

function gate_x(sites::Vector{<:Index}, n::Int)
    return adjointmap_itensor("X", sites, n)
end

function gate_y(sites::Vector{<:Index}, n::Int)
    return adjointmap_itensor("Y", sites, n)
end

function gate_z(sites::Vector{<:Index}, n::Int)
    return adjointmap_itensor("Z", sites, n)
end

function gate_p(sites::Vector{<:Index}, n::Int, ϕ::Real)
    return adjointmap_itensor("Phase", sites, n; ϕ=ϕ)
end

function gate_cp(sites::Vector{<:Index}, control::Int, target::Int, ϕ::Real)
    return adjointmap_itensor("CPhase", sites, control, target; ϕ=ϕ)
end

function gate_s(sites::Vector{<:Index}, n::Int)
    return adjointmap_itensor("S", sites, n)
end

function gate_sdg(sites::Vector{<:Index}, n::Int)  # Adjoint of "S"
    # S = Phase(-π/2), where
    #
    #             ⎛ 1       0      ⎞
    # Phase(ϕ) =  ⎜                ⎟
    #             ⎝ 0  exp(im * ϕ) ⎠
    #
    # so S* = Phase(π/2)* = Phase(-π/2)
    return adjointmap_itensor("Phase", sites, n; ϕ=-pi / 2)
end

function gate_h(sites::Vector{<:Index}, n::Int)
    return adjointmap_itensor("H", sites, n)
end

function gate_t(sites::Vector{<:Index}, n::Int)
    return adjointmap_itensor("T", sites, n)
end

function gate_tdg(sites::Vector{<:Index}, n::Int)  # Adjoint of "T"
    # T = Phase(π/4), so T* = Phase(π/4)* = Phase(-π/4)
    return adjointmap_itensor("Phase", sites, n; ϕ=-pi / 4)
end

function gate_ccx(sites::Vector{<:Index}, control1::Int, control2::Int, target::Int)
    return adjointmap_itensor("Toffoli", sites, control1, control2, target)
end

function gate_c3x(
    sites::Vector{<:Index}, control1::Int, control2::Int, control3::Int, target::Int
)  # CCCNOT
    return adjointmap_itensor("CCCNOT", sites, control1, control2, control3, target)
end

function gate_c4x(
    sites::Vector{<:Index},
    control1::Int,
    control2::Int,
    control3::Int,
    control4::Int,
    target::Int,
)  # CCCCNOT
    return adjointmap_itensor(
        "CCCCNOT", sites, control1, control2, control3, control4, target
    )
end

function gate_rx(sites::Vector{<:Index}, n::Int, θ::Number)
    return adjointmap_itensor("Rx", sites, n; θ=θ)
end

function gate_ry(sites::Vector{<:Index}, n::Int, θ::Number)
    return adjointmap_itensor("Ry", sites, n; θ=θ)
end

function gate_rz(sites::Vector{<:Index}, n::Int, θ::Number)
    return adjointmap_itensor("Rz", sites, n; θ=θ)
end

function gate_cy(sites::Vector{<:Index}, control::Int, target::Int)
    return adjointmap_itensor("CY", sites, control, target)
end

function gate_cz(sites::Vector{<:Index}, control::Int, target::Int)
    return adjointmap_itensor("CZ", sites, control, target)
end

function gate_ch(sites::Vector{<:Index}, control::Int, target::Int)
    return adjointmap_itensor("CH", sites, control, target)
end

function gate_swap(sites::Vector{<:Index}, n1::Int, n2::Int)
    return adjointmap_itensor("Swap", sites, n1, n2)
end

function gate_cswap(sites::Vector{<:Index}, control::Int, target1::Int, target2::Int)
    return adjointmap_itensor("CSwap", sites, control, target1, target2)
end

function gate_crx(sites::Vector{<:Index}, control::Int, target::Int, θ::Number)
    return adjointmap_itensor("CRx", sites, control, target; θ=θ)
end

function gate_cry(sites::Vector{<:Index}, control::Int, target::Int, θ::Number)
    return adjointmap_itensor("CRy", sites, control, target; θ=θ)
end

function gate_crz(sites::Vector{<:Index}, control::Int, target::Int, λ::Number)
    return adjointmap_itensor("CRz", sites, control, target; θ=λ)
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

function gate_cu1(sites::Vector{<:Index}, control::Int, target::Int, λ::Number)
    return adjointmap_itensor("CU1", sites, control, target; λ=λ)
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
    # It gives |0⟩⟨0| ⊗ I₂ + |1⟩⟨1| ⊗ U₁(λ) where
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
    return adjointmap_itensor("CU3", sites, control, target; θ=θ, ϕ=ϕ, λ=λ)
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
