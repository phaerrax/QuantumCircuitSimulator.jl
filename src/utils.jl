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
    qiskitvector(v::MPS)

Return the coefficient vector of `v` in the computational basis of the n-th tensor product
of ``ℂ²`` with itself, with the basis ordered as Qiskit does.
"""
function qiskitvector(v::MPS)
    n = length(v)
    return [coefficient(v, i; qiskit=true) for i in 0:(2^n - 1)]
end
