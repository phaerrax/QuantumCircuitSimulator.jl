import ITensors: op

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
Base.isequal(p::PauliString, q::PauliString) = isequal(p.string, q.string)
Base.reverse(p::PauliString) = PauliString(reverse(p.string))

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
    # Reverse the Pauli string to comply with Qiskit's convention
    rp = reverse(p)
    x = ITensors.OneITensor()
    for (s, i) in zip(string.(pauli_inttochar.(operators(rp))), indices(rp))
        x *= op(sites, s, i)
    end
    return x
end
