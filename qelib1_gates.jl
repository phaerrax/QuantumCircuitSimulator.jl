# OpenQASM --> ITensor gate map
# -----------------------------
# 
# The gate set is taken from the Quantum Experience standard header
using ITensors

# We need some operators that ITensors does not define.
function ITensors.op(::OpName"U", ::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
  # Return a generic single-qbit SU(2) gate, as defined in arXiv:1707.03429v2.
  #
  # U (θ, ϕ, λ) := Rz(ϕ) Ry(θ) Rz(λ) =
  #
  #     ⎛ ℯ^{−i(ϕ+λ)/2} cos(θ/2)   −ℯ^{−i(ϕ−λ)/2} sin(θ/2) ⎞
  #   = ⎜                                                  ⎟
  #     ⎝ ℯ^{i(ϕ−λ)/2} sin(θ/2)    ℯ^{i(ϕ+λ)/2} cos(θ/2)   ⎠
  return [
          ℯ^(−im*(ϕ+λ)/2)*cos(θ/2)   −ℯ^(−im*(ϕ−λ)/2)*sin(θ/2)
          ℯ^(im*(ϕ−λ)/2)*sin(θ/2)     ℯ^(im*(ϕ+λ)/2)*cos(θ/2)  
         ]
end

function ITensors.op(::OpName"CCCCNOT", st::SiteType"Qubit")
    id = op(OpName("id"), st)
    not = op(OpName("X"),st)
    return kron(id, id, id, id, not)
end

function ITensors.op(::OpName"CH", st::SiteType"Qubit")
    id = op(OpName("id"), st)
    h = op(OpName("H"),st)
    return kron(id, h)
end

ITensors.op(::OpName"id", ::SiteType"Qubit") = [
                                                1 0
                                                0 1
                                               ]

function gate_id(sites::Vector{<:Index}, n::Int)
    return ITensors.op("id", sites, n)
end

function gate_u1(sites::Vector{<:Index}, n::Int; λ::Real)
    return ITensors.op("U", sites, n; θ=0, ϕ=0, λ=λ)
end

function gate_u2(sites::Vector{<:Index}, n::Int; ϕ::Real, λ::Real)
    return ITensors.op("U", sites, n; θ=pi/2, ϕ=ϕ, λ=λ)
end

function gate_u3(sites::Vector{<:Index}, n::Int; θ::Real, ϕ::Real, λ::Real)
    return ITensors.op("U", sites, n; θ=θ, ϕ=ϕ, λ=λ)
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
    return gate_u1(sites, n;  λ=pi/2)
end

function gate_sdg(sites::Vector{<:Index}, n::Int)  # Adjoint of "S"
    return gate_u1(sites, n;  λ=-pi/2)
end

function gate_h(sites::Vector{<:Index}, n::Int)
    return ITensors.op("H", sites, n)
end

function gate_t(sites::Vector{<:Index}, n::Int)
    return gate_u1(sites, n;  λ=pi/4)
end

function gate_tdg(sites::Vector{<:Index}, n::Int)
    return gate_u1(sites, n;  λ=-pi/4)
end

function gate_ccx(sites::Vector{<:Index}, control1::Int, control2::Int, target::Int)
    return ITensors.op("Toffoli", sites, control1, control2, target)
end

function gate_c3x(sites::Vector{<:Index}, control1::Int, control2::Int, control3::Int, target::Int)  # CCCNOT
    return ITensors.op("CCCNOT", sites, control1, control2, control3, target)
end

function gate_c4x(sites::Vector{<:Index}, control1::Int, control2::Int, control3::Int, control4::Int, target::Int)  # CCCCNOT
    return ITensors.op("CCCCNOT", sites, control1, control2, control3, control4, target)
end

function gate_rx(sites::Vector{<:Index}, n::Int; θ::Number)
    return ITensors.op("Rx", sites, n; θ=θ)
end

function gate_ry(sites::Vector{<:Index}, n::Int; θ::Number)
    return ITensors.op("Ry", sites, n; θ=θ)
end

function gate_rz(sites::Vector{<:Index}, n::Int; θ::Number)
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

function gate_crx(sites::Vector{<:Index}, control::Int, target::Int; θ::Number)
    return ITensors.op("CRx", sites, control, target; θ=θ)
end

function gate_cry(sites::Vector{<:Index}, control::Int, target::Int; θ::Number)
    return ITensors.op("CRy", sites, control, target; θ=θ)
end

function gate_crz(sites::Vector{<:Index}, control::Int, target::Int; λ::Number)
    #return ITensors.op("CRz", sites, control, target; θ=λ)
    # ITensors already provides a `CRz` gate, but we follow the implementation in
    # the qelib1.inc file anyway.
    return gate_u1(sites, target; λ=λ/2) * gate_cx(sites, control, target) * gate_u1(sites, target; λ=-λ/2) * gate_cx(sites, control, target)
end

function gate_cu1(sites::Vector{<:Index}, control::Int, target::Int; λ::Number)
    return gate_u1(sites, control, λ=λ/2) *gate_cx(sites, target, control) *gate_u1(sites, target, λ=-λ/2) *gate_cx(sites, target, control) *gate_u1(sites, target, λ=λ/2)
end

function gate_cu3(sites::Vector{<:Index}, control::Int, target::Int;θ::Real, ϕ::Real, λ::Real)
    # Implements a controlled U(theta,phi,lambda) gate.
    return gate_u1(sites, target; λ=(λ-ϕ)/2) *gate_cx(sites, control, target) *gate_u3(sites, target; θ=-θ/2,ϕ=0,λ=-(ϕ+λ)/2) *gate_cx(sites, control, target) *gate_u3(sites, target; θ=θ/2,ϕ=ϕ,λ=0)
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
