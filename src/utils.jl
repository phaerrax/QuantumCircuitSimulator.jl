using ITensors

"""
    coefficient(v::MPS, i::Integer; qiskit=false)

Return the `i`-th coefficient of `v` in the computational basis of the n-th tensor
product of ``ℂᵈ`` with itself, with ``d`` being the local dimension. Optionally use
Qiskit's way of ordering the basis.
"""
function coefficient(v::MPS, i::Integer; qiskit=false)
    s = noprime(siteinds(first, v))
    !allequal(dim.(s)) && error("Local dimensions differ between sites.")
    locdim = first(dim.(s))

    n = length(v)
    inds = digits(i; base=locdim, pad=n)  # e.g. digits(3; base=2, pad=3) == [1, 1, 0]
    if qiskit
        # Qiskit orders the computational basis as 00 01 10 11, and also the last digit
        # is the first qbit in the order. An "XY" operator with X on qbit 0 and Y on qbit 1
        # is actually represented as Y ⊗ X.
        reverse!(inds)
    end
    vec = ITensors.OneITensor()
    for j in 1:length(v)
        vec *= (v[j] * state(s[j], inds[j] + 1))
        # Julia indices are 1, 2, ... ------^
    end
    return scalar(vec)
end

"""
    coefficient(m::MPO, i::Integer, j::Integer; qiskit=false)

Return the `(i,j)`-th coefficient of `m` in the computational basis of the n-th tensor
product of ``ℂᵈ`` with itself, with ``d`` being the local dimension. Optionally use
Qiskit's way of ordering the basis.
"""
function coefficient(v::MPO, i::Integer, j::Integer; qiskit=false)
    s = noprime(siteinds(first, v))
    !allequal(dim.(s)) && error("Local dimensions differ between sites.")
    locdim = first(dim.(s))

    n = length(v)
    inds_i = digits(i; base=locdim, pad=n) # e.g. digits(3; base=2, pad=3) == [1, 1, 0]
    inds_j = digits(j; base=locdim, pad=n)
    # For example, with locdim=4 we get
    #   inds_i = [
    #      [0, 0, 0, 0, 0]
    #      [1, 0, 0, 0, 0]
    #      [2, 0, 0, 0, 0]
    #      [3, 0, 0, 0, 0]
    #      [0, 1, 0, 0, 0]
    #      [1, 1, 0, 0, 0]
    #      [2, 1, 0, 0, 0]
    #      …
    #   ]
    # so the most significant digits are at higher indices, such that
    #   n == sum(digits[k]*base^(k-1) for k=1:length(digits))

    if qiskit
        # Our standard order is 00 10 01 11.
        # Qiskit orders the computational basis as 00 01 10 11.
        reverse!(inds_i)
        reverse!(inds_j)
    end
    vec = ITensors.OneITensor()
    for k in 1:length(v)
        vec *= (dag(state(s[k]', inds_i[k] + 1)) * v[k] * state(s[k], inds_j[k] + 1))
        # Julia indices are 1, 2, ... ------^
    end
    return scalar(vec)
end

basis_order_message = """
The basis is ordered with the most significant digit at higher indices; for example
if ``d=2`` and ``n=3`` it is ordered as

```
[0, 0, 0]
[1, 0, 0]
[0, 1, 0]
[1, 1, 0]
[0, 0, 1]
[1, 0, 1]
[0, 1, 1]
[1, 1, 1]
```

i.e. the ``i``-th element of the basis is ``e_{i_1} ⊗ e_{i_2} ⊗ ⋯ ⊗ e_{i_n}`` where
``(i_1, …, i_n)`` are the digits in `digits(i-1, base=d, pad=n)`.

Optionally the "Qiskit order", where the digit sequences are reversed, so for example
`[0, 0, 0]` would be followed by `[0, 0, 1]` in the examples above, can be switched on by
the appropriate keyword argument.
"""

"""
    fullvector(v::MPS; qiskit=false)

Return the coefficient vector of `v` in the computational basis of the ``n``-th tensor
product of ``ℂᵈ`` with itself, ``d`` being the local dimension and ``n`` the number of
sites.
$basis_order_message
"""
function fullvector(v::MPS; qiskit=false)
    dims = dim.(siteinds(first, v))
    allequal(dims) || error("Not all sites of the MPS have the same dimension.")
    d = dims[1]

    n = length(v)
    vec = Vector{ComplexF64}(undef, d^n)
    for i in 1:(d^n)
        vec[i] = coefficient(v, i - 1; qiskit=qiskit)
    end
    return vec
end

"""
    fullmatrix(m::MPO; qiskit=false)

Return the coefficient matrix of `m` in the computational basis of the ``n``-th tensor
product of ``ℂᵈ`` with itself, ``d`` being the local dimension and ``n`` the number of
sites.
$basis_order_message
"""
function fullmatrix(m::MPO; qiskit=false)
    dims = dim.(siteinds(first, m))
    allequal(dims) || error("Not all sites of the MPS have the same dimension.")
    d = dims[1]

    n = length(m)
    mat = Matrix{ComplexF64}(undef, d^n, d^n)
    for i in 1:(d^n), j in 1:(d^n)
        mat[i, j] = coefficient(m, i - 1, j - 1; qiskit=qiskit)
    end
    return mat
end
