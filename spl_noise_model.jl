using OffsetArrays, LinearAlgebra, ITensors

"""
    noise_ptm_coefficients(λ, λλ)

Return the parameters of the noise model in the sparse Pauli-Lindblad model from the
given set of coefficients `λ` and `λλ`, respectively a set of ``N`` 3-element vectors
and a set of ``N`` 3x3 matrices, where ``N`` represents the number of qbits in the
circuit.
"""
function noise_ptm_coefficients(λ, λλ)
    if length(λ) != length(λλ)
        error("Uneven lenghts in parameter arrays")
    end
    n_qbits = length(λ)
    ϰ = [OffsetArray(ones(4), 0:3) for _ in 1:n_qbits]
    ϰϰ = [OffsetArray(ones(4, 4), 0:3, 0:3) for _ in 1:n_qbits]
    for m in 1:n_qbits
        for i in 1:3
            ϰ[m][i] = exp(-2 * sum([λ[m][j] for j in 1:3 if j != i]))
            # sum(A, dims=1) --> sum along columns
            # sum(A, dims=2) --> sum along rows
            #
            # If A = (a b // c d) then
            #   sum(A, dims=1) = (a+c b+d)
            #   sum(A, dims=2) = (a+b // c+d)
            ϰϰ[m][0, i] = exp(-2 * (sum(λλ[m]) - sum(λλ[m]; dims=1)[i]))
            ϰϰ[m][i, 0] = exp(-2 * (sum(λλ[m]) - sum(λλ[m]; dims=2)[i]))
            for j in 1:3
                ϰϰ[m][i, j] = exp(
                    -2 * (sum(λλ[m]; dims=2)[i] + sum(λλ[m]; dims=1)[j] - 2 * λλ[m][i, j])
                )
            end
        end
    end
    return ϰ, ϰϰ
end

# Naive MPO construction
# ----------------------
"""
    noise_inverse_mpo(sites::Vector{<:Index}, λ, λλ)

Return the MPO form of the inverse of the noise map generated by the parameter sets
`λ` and `λλ`, as explained in [1].

# References
[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function noise_inverse_mpo(s::Vector{<:Index}, λ, λλ)
    if !(length(λ) == length(λλ) == length(s))
        error("Uneven lenghts in parameter arrays")
    end
    n_qbits = length(λ)

    ϰ, ϰϰ = noise_ptm_coefficients(λ, λλ)

    # Single-qbit gates: just put the ϰ[m] vector on the diagonal.
    q1g_ptm_matrix = [diagm(0 => ϰ[m]) for m in 1:n_qbits]
    q1g = [ITensor(q1g_ptm_matrix[m], s[m]', dag(s[m])) for m in 1:n_qbits]

    # Two-qbit gates: unravel the ϰϰ[m] matrix by rows, i.e. obtain the vector
    #   (ϰϰ[m][0,0], ϰϰ[m][0,1], ϰϰ[m][0,2], ϰϰ[m][0,3], ϰϰ[m][1,0], ...)
    # and put on the diagonal.
    q2g_ptm_matrix = [diagm(0 => [ϰϰ[m]'...]) for m in 1:n_qbits]
    q2g = [
        ITensor(q2g_ptm_matrix[m], s[m + 1]', s[m]', dag(s[m + 1]), dag(s[m])) for
        m in 1:(n_qbits - 1)
    ]

    invN = MPO(q1g)
    # This is just the single-qbit gates. Now we apply all the odd gates, and then
    # the even ones.
    for m in 1:2:(n_qbits - 1)
        invN = apply(q2g[m], invN)
    end
    for m in 2:2:(n_qbits - 1)
        invN = apply(q2g[m], invN)
    end

    # We know from the cited paper that the expected bond dimension of the MPO is 4,
    # so we happily truncate it without fear of any consequence (right??).
    return truncate(invN; maxdim=4)
end