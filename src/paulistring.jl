using ITensors.SiteTypes: _sitetypes, commontags
using LindbladVectorizedTensors: vop

import ITensorMPS: MPS

isvalidpauliint(i::Integer) = 0 <= i <= 3

function pauli_chartoint(c::Char)
    if c == 'I'
        return 0
    elseif c == 'X'
        return 1
    elseif c == 'Y'
        return 2
    elseif c == 'Z'
        return 3
    else
        error("Invalid Pauli char")
    end
end

function pauli_inttochar(c::Integer)
    if c == 0
        return 'I'
    elseif c == 1
        return 'X'
    elseif c == 2
        return 'Y'
    elseif c == 3
        return 'Z'
    else
        error("Invalid Pauli integer")
    end
end

struct PauliString
    string::Vector{Int8}
    function PauliString(ints::Vector{<:Integer})
        return all(isvalidpauliint, ints) ? new(ints) : error("Invalid Pauli integer")
    end
end

"""
    PauliString(L::Integer, t::Tuple{Char, Integer}...)

Return a PauliString of length `L` with the given Pauli operators.

# Example

```julia-repl
julia> p = PauliString(8, ('X', 3), ('Y', 4), ('Z', 6))
PauliString(Int8[0, 0, 1, 2, 0, 3, 0, 0])

julia> show(p)
"IIXYIZII"
```
"""
function PauliString(L::Integer, ts::Tuple{Char,Integer}...)
    str = repeat([0], L)
    for t in ts
        str[t[2]] = pauli_chartoint(t[1])
    end
    return PauliString(str)
end

function PauliString(str::AbstractString)
    ints = Vector{Int8}(undef, length(str))
    for i in eachindex(str)
        ints[i] = pauli_chartoint(str[i])
    end
    return PauliString(ints)
end

# Overload common methods for PauliStrings

Base.string(p::PauliString) = join(pauli_inttochar.(p.string))
Base.length(p::PauliString) = length(p.string)
Base.show(p::PauliString) = show(string(p))
Base.print(p::PauliString) = print(string(p))
Base.getindex(p::PauliString, i::Integer) = getindex(p.string, i)
function Base.setindex!(p::PauliString, c::Integer, i::Integer)
    if isvalidpaulichar(c)
        return setindex!(p.string, c, i)
    else
        error("Invalid Pauli integer")
    end
end
Base.setindex!(p::PauliString, c::Char, i::Integer) = setindex!(p, pauli_chartoint(c), i)

Base.isless(p::PauliString, q::PauliString) = isless(p.string, q.string)

# - Equality testing -
# Since our Pauli strings are just arrays of integers, we fall back to the existing
# implementations for that type.
Base.isequal(p::PauliString, q::PauliString) = isequal(p.string, q.string)
# We need to implement the `hash` function as well, otherwise functions such as `unique`
# will not work.
Base.:(==)(p::PauliString, q::PauliString) = p.string == q.string
Base.hash(p::PauliString, h::UInt) = hash(p.string, h)

Base.reverse(p::PauliString) = PauliString(reverse(p.string))

"""
    istrivial(p::PauliString)

Return `true` if all factors of `p` are the identity.
"""
istrivial(p::PauliString) = all(==(pauli_chartoint('I')), p.string)

"""
    crop(p::PauliString, range)

Return a subset of the Pauli string `p` containing only the factors whose index is in
`range` (which must be sorted).
"""
function crop(p::PauliString, range)
    # Sanity checks
    !issorted(range) && error("Range $range is not sorted")
    !issubset(range, eachindex(p.string)) && error("Range $range out of bounds")

    return PauliString(p.string[range])
end

"""
    indices(p::PauliString)

Return the site indices of non-trivial factors in `p`.
"""
function indices(p::PauliString)
    return findall(!=(0), p.string)
end

"""
    order(p::PauliString)

Return the number of non-trivial factors in the Pauli string `p`.

# Example

```julia-repl
julia> p = PauliString(8, ('X', 3), ('Y', 4), ('Z', 6))
PauliString(Int8[0, 0, 1, 2, 0, 3, 0, 0])

julia> order(p)
3
```
"""
order(p::PauliString) = length(indices(p))

"""
    operators(p::PauliString)

Return a list containing the non-trivial factors in `p` (ordered).
"""
function operators(p::PauliString)
    return filter(!=(0), p.string)
end

"""
    ITensors.op(p::PauliString)

Return an ITensor corresponding the Pauli string operator.
"""
function ITensors.op(sites::Vector{<:Index}, p::PauliString)
    length(p) != length(sites) && "Lengths of Pauli string and Index vector differ."
    x = ITensors.OneITensor()
    for (s, i) in zip(string.(pauli_inttochar.(operators(p))), indices(p))
        x *= op(sites, s, i)
    end
    return x
end

ITensorMPS.MPS(::SiteType, sites::Vector{<:Index}, p::PauliString) = nothing
ITensorMPS.MPO(::SiteType, sites::Vector{<:Index}, p::PauliString) = nothing

"""
    MPS(sites::Vector{<:Index}, p::PauliString)

Return an MPS corresponding to the Pauli string `p` on the site indices `sites`.
"""
function ITensorMPS.MPS(sites::Vector{<:Index}, p::PauliString)
    length(p) != length(sites) && "Lengths of Pauli string and Index vector differ."
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    for st in common_stypes
        res = ITensorMPS.MPS(st, sites, p)
        !isnothing(res) && return res
    end

    return throw(
        ArgumentError(
            "Overload of \"MPS\" function not found for gate name \"$name\" and Index " *
            "tags: $(tags.(sites)).",
        ),
    )
end

"""
    MPO(sites::Vector{<:Index}, p::PauliString)

Return an MPO corresponding to the Pauli string `p` on the site indices `sites`.
"""
function ITensorMPS.MPO(sites::Vector{<:Index}, p::PauliString)
    length(p) != length(sites) && "Lengths of Pauli string and Index vector differ."
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    for st in common_stypes
        res = ITensorMPS.MPO(st, sites, p)
        !isnothing(res) && return res
    end

    return throw(
        ArgumentError(
            "Overload of \"MPO\" function not found for gate name \"$name\" and Index " *
            "tags: $(tags.(sites)).",
        ),
    )
end

function ITensorMPS.MPS(::SiteType"vQubit", sites::Vector{<:Index}, p::PauliString)
    statenames = string.(pauli_inttochar.(p.string))
    statenames = replace.(statenames, "I" => "Id")
    vstatenames = ["v$sn" for sn in statenames]
    return MPS(sites, vstatenames)
end

function ITensorMPS.MPO(::SiteType"Qubit", sites::Vector{<:Index}, p::PauliString)
    opnames = string.(pauli_inttochar.(p.string))
    opnames = replace.(opnames, "I" => "Id")
    return MPO(ComplexF64, sites, opnames)
end

"""
    samplepaulistrings(v::MPS, nsamples::Integer)

Sample `nsamples` Pauli strings from `v` and compute their overlap with the MPS.
Return a pair `ps, overlaps` where `overlaps[k]` is the coefficient of the `ps[k]`
component of `v`: this means that if we write `v` as a linear combination of Pauli
strings ``v = ∑ₖ cₖσₖ`` then `overlaps[k]` is the coefficient ``cₖ``.

Note that the MPS of a Pauli string is not normalized in the Hilbert-Schmidt inner
product ``⟨A,B⟩ = tr(A† B)``: the norm of a Pauli string of length `N` is ``2^(N/2)``.
"""
function samplepaulistrings(v::MPS, nsamples::Integer)
    if any(t -> !(SiteType("vQubit") in t), _sitetypes.(siteinds(x)))
        error("samplepaulistrings works for vQubit site types only.")
    end

    # The sampling algorithm requires a normalized MPS orthogonalized on the first site.
    vn = orthogonalize(v, 1)
    vn /= norm(vn)
    sites = siteinds(v)

    ps = Vector{PauliString}(undef, nsamples)
    for i in 1:nsamples
        ps[i] = PauliString(sample(vn) .- 1)
        # sample(vn) gives us a vector of elements from {1, 2, 3, 4}; they are indices
        # referring to the PTM basis {I/sqrt(2), X/sqrt(2), Y/sqrt(2), Z/sqrt(2)}.
        # PauliString objects work with {0, 1, 2, 3} (respectively) instead so we need to
        # decrease by one.
    end

    # We compute the overlap of each Pauli string with the (normalized) MPS v.
    # Note that the MPS created from a PauliString object is not normalized, since it uses
    # the {I, X, Y, Z} matrices instead of their rescaled counterparts in the PTM basis.
    # For this reason, each MPS(sites, p) has norm 2^(length(sites)/2), and the overlap
    # must be normalized accordingly.
    # The orthonormal PTM basis is (N := length(sites))
    #   eₖ = 2^(-N/2) σₖ,
    # so we can extract the coefficients with the inner product as follows:
    #   v = ∑ₖ cₖσₖ
    #   ⟨eⱼ, v⟩ = ∑ₖ cₖ⟨eⱼ,σₖ⟩ = 2^(N/2) ∑ₖ cₖ⟨eⱼ,eₖ⟩ =  2^(N/2) cⱼ
    # therefore cₖ = 2^(-N/2) ⟨eₖ, v⟩ = 2^(-N) ⟨σₖ, v⟩.
    normf = 2^length(sites)
    overlaps = [dot(MPS(sites, p), v) / normf for p in ps]  # FIXME cache the contractions
    return ps, overlaps
end
