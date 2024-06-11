using ITensors.SiteTypes: _sitetypes, commontags

export gate

# We cannot directly use ITensors' own "op" here for our gates since we would like to define
# gates with the same name acting on Qubit and vQubit site types alike, the only difference
# being in the site type.
# For example, we would like to have a gate "U" acting as the operator "U" on qbits and as
# "U* ⋅ U" on vectorized qbits. We can't define the second one as "U" for vQubits in a way
# that is compatible with the automatic rules built in LindbladVectorizedTensors to build
# multiplication operators.

# Default implementations of gate
# -------------------------------
gate(::GateName; kwargs...) = nothing
gate(::GateName, ::SiteType; kwargs...) = nothing
gate(::GateName, ::SiteType, ::Index...; kwargs...) = nothing
function gate(
    ::GateName, ::SiteType, ::SiteType, sitetypes_inds::Union{SiteType,Index}...; kwargs...
)
    return nothing
end

# Main definition of gate
# -----------------------
"""
    gate(name::String, s::Index...; kwargs...)

Return an ITensor corresponding to the gate named `name` for the Index `s`.
The operator is constructed by calling an overload of the `gate` method which takes a
`SiteType` argument that corresponds to one of the tags of the Index `s` and an
`GateName"name"` argument that corresponds to the input gate name.
If the gate requires numerical parameters (a.k.a. "classical arguments") too, they must be
provided within a keyword argument named `cargs`.

# Example

```julia
s = siteinds("Qubit", 2)
z = gate("z", s[1])
u3 = gate("u3", s[1]; cargs=(0, pi/4, 2))
cp = gate("cp", s[1], s[2]; cargs=(pi/2))
```
"""
function gate(name::AbstractString, s::Index...; kwargs...)
    # This is a bare-bones version of ITensors.op, stripped of many features which we don't
    # need here. We also removed all variants which ITensors keeps for compatibility with
    # previous versions.
    name = strip(name)
    gn = GateName(name)
    commontags_s = commontags(s...)
    common_stypes = _sitetypes(commontags_s)

    # (We skip all the algebra machinery ITensors uses in its `op` function, which we
    # will not support here. A gate name is just a gate name, no operations allowed in the
    # input strings.)

    # 1) Try calling a function of the form:
    #    gate(::GateName, ::SiteType, ::Index...; kwargs...)
    # which already returns an ITensor
    for st in common_stypes
        res = gate(gn, st, s...; kwargs...)
        !isnothing(res) && return res
    end

    # 2) otherwise try calling a function of the form:
    #    gate(::GateName, ::SiteType; kwargs...)
    # which returns a Julia matrix
    for st in common_stypes
        g_mat = gate(gn, st; kwargs...)
        if !isnothing(g_mat)
            rs = reverse(s)
            return itensor(g_mat, prime.(rs)..., ITensors.dag.(rs)...)
        end
    end

    # (Here ITensors.op starts looking for overloads for common tags found, in case the caller
    # is trying to define an operator with mixed site types. We don't need it here.)
    return throw(
        ArgumentError(
            "Overload of \"gate\" function not found for gate name \"$name\" and Index tags: $(tags.(s)).",
        ),
    )
end

gate(name::AbstractString; kwargs...) = error("Must input indices when creating a `gate`.")

# On-the-fly gate construction
# ----------------------------
# Useful for creating a new gate which is not a predefined one.
"""
    gate(X::AbstractArray, s::Index...)
    gate(M::Matrix, s::Index...)

Given a matrix `M` and a set of indices `s`, `t`, ... return a gate ITensor with matrix
elements given by `M` and indices `s, s', t, t'`.

# Examples

```julia
julia> s = siteind("Qubit")
(dim=2|id=575|"Qubit,Site")

julia> g = gate([1/2 0; 0 -1/2], s)
ITensor ord=2 (dim=2|id=575|"Qubit,Site")' (dim=2|id=575|"Qubit,Site")
NDTensors.Dense{Float64, Vector{Float64}}

julia> @show g
g = ITensor ord=2
Dim 1: (dim=2|id=575|"Qubit,Site")'
Dim 2: (dim=2|id=575|"Qubit,Site")
NDTensors.Dense{Float64, Vector{Float64}}
 2×2
 0.5   0.0
 0.0  -0.5
ITensor ord=2 (dim=2|id=575|"Qubit,Site")' (dim=2|id=575|"Qubit,Site")
NDTensors.Dense{Float64, Vector{Float64}}
```
"""
gate(X::AbstractArray, s::Index...) = itensor(X, prime.([s...]), dag.([s...]))

# Splat vectors of indices
gate(gatename, s::Vector{<:Index}; kwargs...) = gate(gatename, s...; kwargs...)
gate(s::Vector{<:Index}, gatename; kwargs...) = gate(gatename, s...; kwargs...)

# To ease calling of other gate overloads, allow passing a string as the gate name
function gate(gatename::AbstractString, t::SiteType; kwargs...)
    return gate(GateName(gatename), t; kwargs...)
end

# Gate constructors with name, a list of sites and some integers that point to which sites
# the gate acts on.
"""
    gate(gatename, s::Vector{<:Index}, ns::NTuple{N,Integer}; kwargs...)

Return an ITensor corresponding to the gate named `gatename` on sites `s[n]` for each `n`
in the tuple `ns`.

# Example

```julia
s = siteinds("Qubit", 4)

g2 = gate("x", s, 2)
g13 = gate("cnot", s, (1, 3))
```
"""
function gate(gatename, s::Vector{<:Index}, ns::NTuple{N,Integer}; kwargs...) where {N}
    return gate(gatename, ntuple(n -> s[ns[n]], Val(N))...; kwargs...)
    # ntuple(f, ::Val{N}) creates a tuple of length N, computing each element as f(i)
end

"""
    gate(gatename, s::Vector{<:Index}, ns::Vararg{Integer}; kwargs...)

Return an ITensor corresponding to the gate named `gatename` on sites `s[n]` for each `n`
in `ns`.

# Example

```julia
s = siteinds("Qubit", 4)

g2 = gate("x", s, 2)
g13 = gate("cnot", s, 1, 3)
```
"""
function gate(gatename, s::Vector{<:Index}, ns::Vararg{Integer}; kwargs...)
    return gate(gatename, s, ns; kwargs...)
end

# The following are all slight variations on the two methods above.

function gate(s::Vector{<:Index}, gatename, ns::Tuple{Vararg{Integer}}; kwargs...)
    # Vararg{T} corresponds to zero or more elements of type T.
    return gate(gatename, s, ns...; kwargs...)
end

function gate(s::Vector{<:Index}, gatename, ns::Integer...; kwargs...)
    return gate(gatename, s, ns; kwargs...)
end

function gate(s::Vector{<:Index}, gatename, ns::Tuple{Vararg{Integer}}, kwargs::NamedTuple)
    return gate(gatename, s, ns; kwargs...)
end

function gate(s::Vector{<:Index}, gatename, ns::Integer, kwargs::NamedTuple)
    return gate(gatename, s, (ns,); kwargs...)
end

gate(s::Vector{<:Index}, o::Tuple) = gate(s, o...)

gate(o::Tuple, s::Vector{<:Index}) = gate(s, o...)
