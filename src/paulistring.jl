using ITensors.SiteTypes: _sitetypes, commontags
using LindbladVectorizedTensors: vop

import ITensorMPS: MPS, MPO

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
    # String sanity check:
    # 1) all numbers shall be between 1 and L (inclusive)
    if !all(1 .<= last.(ts) .<= L)
        error("Some Pauli operators do not fit in the given length.")
    end
    # 2) no char shall be repeated
    if !allunique(last.(ts))
        error("Repeated Pauli operators.")
    end

    str = repeat([0], L)
    for t in ts
        str[t[2]] = pauli_chartoint(t[1])
    end
    return PauliString(str)
end

"""
    PauliString(L::Integer, str::AbstractString)

Return a PauliString of length `L`, with Pauli operators specified in `str` as a sequence
of letters and numbers.

# Example

```julia-repl
julia> p = PauliString(8, "X3Y4Z6")
PauliString(Int8[0, 0, 1, 2, 0, 3, 0, 0])

julia> show(p)
"IIXYIZII"
```
"""
function PauliString(L::Integer, str::AbstractString)
    idx = first.(findall(r"[XYZI]", str))
    factors = [
        [str[idx[i]:(idx[i + 1] - 1)] for i in 1:(length(idx) - 1)]
        str[idx[end]:end]
    ]
    @assert *(factors...) == str
    return PauliString(L, [(first(f), parse(Int, f[2:end])) for f in factors]...)
end

function PauliString(str::AbstractString)
    ints = Vector{Int8}(undef, length(str))
    for i in eachindex(str)
        ints[i] = pauli_chartoint(str[i])
    end
    return PauliString(ints)
end

# Overload common methods for PauliStrings

Base.length(p::PauliString) = length(p.string)

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

# Pretty-printing methods

"""
    compactstring(p::PauliString)

Return the Pauli string in compact form, i.e. written as a sequence of its non-trivial
operators and their position.

# Example

```julia-repl
julia> p = PauliString(8, "X3Y4Z6")
PauliString(Int8[0, 0, 1, 2, 0, 3, 0, 0])

julia> show(p)
"IIXYIZII"

julia> compactstring(p)
"X3Y4Z6"
```
"""
function compactstring(p::PauliString)
    return join(Iterators.flatten(zip(pauli_inttochar.(operators(p)), indices(p))))
end

Base.string(p::PauliString) = join(pauli_inttochar.(p.string))
function Base.show(io::IO, p::PauliString)
    compact = get(io, :compact, false)
    return show(io, compact ? compactstring(p) : string(p))
end

# ITensor interoperability

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

MPS(::SiteType, sites::Vector{<:Index}, p::PauliString) = nothing
MPO(::SiteType, sites::Vector{<:Index}, p::PauliString) = nothing

"""
    MPS(sites::Vector{<:Index}, p::PauliString)

Return an MPS corresponding to the Pauli string `p` on the site indices `sites`.
"""
function MPS(sites::Vector{<:Index}, p::PauliString)
    length(p) != length(sites) && "Lengths of Pauli string and Index vector differ."
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    for st in common_stypes
        res = MPS(st, sites, p)
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
function MPO(sites::Vector{<:Index}, p::PauliString)
    length(p) != length(sites) && "Lengths of Pauli string and Index vector differ."
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    for st in common_stypes
        res = MPO(st, sites, p)
        !isnothing(res) && return res
    end

    return throw(
        ArgumentError(
            "Overload of \"MPO\" function not found for gate name \"$name\" and Index " *
            "tags: $(tags.(sites)).",
        ),
    )
end

function MPS(::SiteType"vQubit", sites::Vector{<:Index}, p::PauliString)
    statenames = string.(pauli_inttochar.(p.string))
    statenames = replace.(statenames, "I" => "Id")
    vstatenames = ["v$sn" for sn in statenames]
    return MPS(sites, vstatenames)
end

function MPO(::SiteType"Qubit", sites::Vector{<:Index}, p::PauliString)
    opnames = string.(pauli_inttochar.(p.string))
    opnames = replace.(opnames, "I" => "Id")
    return MPO(ComplexF64, sites, opnames)
end

@memoize function _contract(x::ITensor, y::ITensor)
    # This tiny little function is here just so that we can memoize the calculations.
    #
    # Note: @memoize by default compares arguments by using the hash function `objectid`.
    # If a primitive type (?) is supplied directly (without first binding it to a variable)
    # then the object id is the same; with composite types this might not happen.
    # For example, if we have a memoized function `f` of a single variable, then repeated
    # calls to f(4) will use the cached value, but f([1, 2, 3]) won't, since [1, 2, 3] is
    # treated as a new variable each time. We should keep that in mind.
    # Anyway, the `shotmeasurement` function is called with the same `shotduals` and
    # `observable` objects each time (each one bound to a variable) so it should be fine.

    # Some stats from @btime:
    #   ⋅ ps = vector of 50 Pauli strings
    #   ⋅ s = siteinds("vQubit", 50)
    #   ⋅ v = randomMPS(s; linkdims=10)
    #
    # @btime foreach(p -> dot(MPS(s, p), v), ps)
    #   216.343 ms (970371 allocations: 233.83 MiB)
    #
    # @btime foreach(p -> QuantumCircuitSimulator.contractpauli(p, v), ps)
    #   30.246 ms (209801 allocations: 36.17 MiB)
    return x * y
end
