using ITensors

"""
    coefficient(v::MPS, i::Integer; qiskit=true)

Return the `i`-th coefficient of `v` in the computational basis of the n-th tensor product
of ``ℂ²`` with itself, with the option to choose if the basis should be ordered as Qiskit
does or the reverse ordering.
"""
function coefficient(v::MPS, i::Integer; qiskit=true)
    s = noprime(siteinds(first, v))
    n = length(v)
    inds = digits(i; base=2, pad=n)
    # for example digits(3; base=2, pad=3) == [1, 1, 0]
    qiskit || reverse!(inds)
    vec = ITensors.OneITensor()
    for j in 1:length(v)
        vec *= (v[j] * state(s[j], inds[j] + 1))
        #      Julia indices are 1, 2 ------^
    end
    return scalar(vec)
end

"""
    coefficient(m::MPO, i::Integer, j::Integer; qiskit=true)

Return the `(i,j)`-th coefficient of `m` in the computational basis of the n-th tensor
product of ``ℂ²`` with itself, with the option to choose if the basis should be ordered
as Qiskit does or the reverse ordering.
"""
function coefficient(v::MPO, i::Integer, j::Integer; qiskit=true)
    s = noprime(siteinds(first, v))
    n = length(v)
    inds_i = digits(i; base=2, pad=n)
    inds_j = digits(j; base=2, pad=n)
    # for example digits(3; base=2, pad=3) == [1, 1, 0]
    qiskit || reverse!(inds)
    vec = ITensors.OneITensor()
    for k in 1:length(v)
        vec *= (dag(state(s[k]', inds_i[k] + 1)) * v[k] * state(s[k], inds_j[k] + 1))
        #   Julia indices are 1, 2 ------^
    end
    return scalar(vec)
end

"""
    qiskitvector(v::MPS)

Return the coefficient vector of `v` in the computational basis of the n-th tensor product
of ``ℂ²`` with itself, with the basis ordered as Qiskit does.
"""
function qiskitvector(v::MPS)
    dims = dim.(siteinds(first, v))
    allequal(dims) || error("Not all sites of the MPS have the same dimension.")
    d = dims[1]

    n = length(v)
    vec = Vector{ComplexF64}(undef, d^n)
    for i in 1:(d^n)
        vec[i] = coefficient(v, i - 1; qiskit=true)
    end
    return vec
end

"""
    qiskitmatrix(m::MPO)

Return the coefficient matrix of `m` in the computational basis of the n-th tensor product
of ``ℂ²`` with itself, with the basis ordered as Qiskit does.
"""
function qiskitmatrix(m::MPO)
    dims = dim.(siteinds(first, m))
    allequal(dims) || error("Not all sites of the MPS have the same dimension.")
    d = dims[1]

    n = length(m)
    mat = Matrix{ComplexF64}(undef, d^n, d^n)
    for i in 1:(d^n), j in 1:(d^n)
        mat[i, j] = coefficient(m, i - 1, j - 1; qiskit=true)
    end
    return mat
end
