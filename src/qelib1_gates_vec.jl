# OpenQASM 3.0 --> ITensor gate map (with vectorized operators)
# -------------------------------------------------------------
# 
# The gate set is taken from the Quantum Experience standard header, and adheres to the
# OpenQASM 3.0 specification.

using LindbladVectorizedTensors: adjointmap_itensor

function gate(::GateName"id", ::SiteType"vQubit", s::Index)
    return ITensors.op("Id", s)
end

function gate(::GateName"u1", ::SiteType"vQubit", s::Index, λ::Real)
    return adjointmap_itensor("U", s; θ=0, ϕ=0, λ=λ)
end

function gate(::GateName"u2", ::SiteType"vQubit", s::Index, ϕ::Real, λ::Real)
    return adjointmap_itensor("U", s; θ=pi / 2, ϕ=ϕ, λ=λ)
end

function gate(::GateName"u3", ::SiteType"vQubit", s::Index, θ::Real, ϕ::Real, λ::Real)
    return adjointmap_itensor("U", s; θ=θ, ϕ=ϕ, λ=λ)
end

function gate(::GateName"u", st::SiteType"vQubit", s::Index, θ::Real, ϕ::Real, λ::Real)
    return gate(GateName("u3"), st, s, θ, ϕ, λ)
end

function gate(::GateName"cx", ::SiteType"vQubit", control::Index, target::Index)
    return adjointmap_itensor("CX", control, target)
end

function gate(::GateName"x", ::SiteType"vQubit", s::Index)
    return adjointmap_itensor("X", s)
end

function gate(::GateName"y", ::SiteType"vQubit", s::Index)
    return adjointmap_itensor("Y", s)
end

function gate(::GateName"z", ::SiteType"vQubit", s::Index)
    return adjointmap_itensor("Z", s)
end

function gate(::GateName"p", ::SiteType"vQubit", s::Index, ϕ::Real)
    return adjointmap_itensor("Phase", s; ϕ=ϕ)
end

function gate(::GateName"cp", ::SiteType"vQubit", control::Index, target::Index, ϕ::Real)
    return adjointmap_itensor("CPhase", control, target; ϕ=ϕ)
end

function gate(::GateName"s", ::SiteType"vQubit", s::Index)
    return adjointmap_itensor("S", s)
end

function gate(::GateName"sdg", ::SiteType"vQubit", s::Index)  # Adjoint of "S"
    # S = Phase(-π/2), where
    #
    #             ⎛ 1       0      ⎞
    # Phase(ϕ) =  ⎜                ⎟
    #             ⎝ 0  exp(im * ϕ) ⎠
    #
    # so S* = Phase(π/2)* = Phase(-π/2)
    return adjointmap_itensor("Phase", s; ϕ=-pi / 2)
end

function gate(::GateName"h", ::SiteType"vQubit", s::Index)
    return adjointmap_itensor("H", s)
end

function gate(::GateName"t", ::SiteType"vQubit", s::Index)
    return adjointmap_itensor("T", s)
end

function gate(::GateName"tdg", ::SiteType"vQubit", s::Index)  # Adjoint of "T"
    # T = Phase(π/4), so T* = Phase(π/4)* = Phase(-π/4)
    return adjointmap_itensor("Phase", s; ϕ=-pi / 4)
end

function gate(
    ::GateName"ccx", ::SiteType"vQubit", control1::Index, control2::Index, target::Index
)
    return adjointmap_itensor("Toffoli", control1, control2, target)
end

function gate(
    ::GateName"c3x",
    ::SiteType"vQubit",
    control1::Index,
    control2::Index,
    control3::Index,
    target::Index,
)
    return adjointmap_itensor("CCCNOT", control1, control2, control3, target)
end

function gate(
    ::GateName"c4x",
    ::SiteType"vQubit",
    control1::Index,
    control2::Index,
    control3::Index,
    control4::Index,
    target::Index,
)
    return adjointmap_itensor("CCCCNOT", control1, control2, control3, control4, target)
end

function gate(::GateName"rx", ::SiteType"vQubit", s::Index, θ::Real)
    return adjointmap_itensor("Rx", s; θ=θ)
end

function gate(::GateName"ry", ::SiteType"vQubit", s::Index, θ::Real)
    return adjointmap_itensor("Ry", s; θ=θ)
end

function gate(::GateName"rz", ::SiteType"vQubit", s::Index, θ::Real)
    return adjointmap_itensor("Rz", s; θ=θ)
end

function gate(::GateName"cy", ::SiteType"vQubit", control::Index, target::Index)
    return adjointmap_itensor("CY", control, target)
end

function gate(::GateName"cz", ::SiteType"vQubit", control::Index, target::Index)
    return adjointmap_itensor("CZ", control, target)
end

function gate(::GateName"ch", ::SiteType"vQubit", control::Index, target::Index)
    return adjointmap_itensor("CH", control, target)
end

function gate(::GateName"swap", ::SiteType"vQubit", s1::Index, s2::Index)
    return adjointmap_itensor("Swap", s1, s2)
end

function gate(
    ::GateName"cswap", ::SiteType"vQubit", control::Index, target1::Index, target2::Index
)
    return adjointmap_itensor("CSwap", control, target1, target2)
end

function gate(::GateName"crx", ::SiteType"vQubit", control::Index, target::Index, θ::Real)
    return adjointmap_itensor("CRx", control, target; θ=θ)
end

function gate(::GateName"cry", ::SiteType"vQubit", control::Index, target::Index, θ::Real)
    return adjointmap_itensor("CRy", control, target; θ=θ)
end

function gate(::GateName"crz", ::SiteType"vQubit", control::Index, target::Index, λ::Real)
    return adjointmap_itensor("CRz", control, target; θ=λ)
end

function gate(::GateName"cu1", ::SiteType"vQubit", control::Index, target::Index, λ::Real)
    return adjointmap_itensor("CU1", control, target; λ=λ)
end

function gate(
    ::GateName"cu3",
    ::SiteType"vQubit",
    control::Index,
    target::Index,
    θ::Real,
    ϕ::Real,
    λ::Real,
)
    return adjointmap_itensor("CU3", control, target; θ=θ, ϕ=ϕ, λ=λ)
end

#function gate(::GateName"c3sqrtx", ::SiteType"vQubit", s::Index)
#    return ITensors.op("", s)
#end
#function gate(::GateName"rxx", ::SiteType"vQubit", s::Index)
#    return ITensors.op("", s)
#end
#
#function gate(::GateName"rzz", ::SiteType"vQubit", s::Index)
#    return ITensors.op("", s)
#end
#
#function gate(::GateName"rccx", ::SiteType"vQubit", s::Index)
#    return ITensors.op("", s)
#end
#
#function gate(::GateName"rc3x", ::SiteType"vQubit", s::Index)
#    return ITensors.op("", s)
#end
