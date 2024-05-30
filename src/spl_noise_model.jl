struct SPLNoiseModel
    parameters::Dict{PauliString,Real}
    function SPLNoiseModel(dict::Dict{PauliString,<:Real})
        # Check that all strings have the same length
        samelength = allequal(length.(keys(dict)))
        # and that their order is 1 or 2
        correctorders = all(n -> (n == 1 || n == 2), [order(p.first) for p in dict])

        if !samelength
            error("Pauli strings don't have the same length")
        elseif !correctorders
            error("Some Pauli strings don't have order 1 or 2")
        else
            new(dict)
        end
    end
end

Base.isempty(model::SPLNoiseModel) = isempty(model.parameters)

function nqbits(model::SPLNoiseModel)
    return isempty(model) ? 0 : length(first(keys(model.parameters)))
end

function SPLNoiseModel(file::AbstractString)
    dict = JSON.parsefile(file)
    return SPLNoiseModel(Dict(PauliString(str) => coeff for (str, coeff) in dict))
end

function noise_ptm_generators(model::SPLNoiseModel)
    N = nqbits(model)

    vec = [zeros(Float64, 3) for _ in 1:N]
    mat = [zeros(Float64, 3, 3) for _ in 1:(N - 1)]
    for (k, v) in model.parameters
        if order(k) == 1
            vec[first(indices(k))][operators(k)...] = v
        elseif order(k) == 2
            mat[first(indices(k))][operators(k)...] = v
        else
            # This shouldn't happens since the inner constructor of SPLNoiseModel enforces
            # this condition... but just in case...
            error("Pauli string $k has order ≠ 1 or 2")
        end
    end

    return vec, mat
end

"""
    noise_ptm_coefficients(ptm_generator_vec, ptm_generator_mat)

Return the parameters of the noise model in the sparse Pauli-Lindblad model from the
given set of coefficients `ptm_generator_vec` and `ptm_generator_mat`, respectively a set
of ``N`` 3-element vectors and a set of ``N`` 3x3 matrices, where ``N`` represents the
number of qbits in the circuit.
"""
function noise_ptm_coefficients(ptm_generator_vec, ptm_generator_mat)
    if length(ptm_generator_vec) == length(ptm_generator_mat) + 1
        nqbits = length(ptm_generator_vec)
        ptm_parameter_vec = [OffsetArray(ones(4), 0:3) for _ in 1:nqbits]
        ptm_parameter_mat = [OffsetArray(ones(4, 4), 0:3, 0:3) for _ in 1:(nqbits - 1)]
        for m in 1:nqbits
            for i in 1:3
                ptm_parameter_vec[m][i] = exp(
                    -2 * sum([ptm_generator_vec[m][j] for j in 1:3 if j != i])
                )
            end
        end
        for m in 1:(nqbits - 1)
            for i in 1:3
                # sum(A, dims=1) --> sum along columns
                # sum(A, dims=2) --> sum along rows
                #
                # If A = (a b // c d) then
                #   sum(A, dims=1) = (a+c b+d)
                #   sum(A, dims=2) = (a+b // c+d)
                ptm_parameter_mat[m][0, i] = exp(
                    -2 * (sum(ptm_generator_mat[m]) - sum(ptm_generator_mat[m]; dims=1)[i])
                )
                ptm_parameter_mat[m][i, 0] = exp(
                    -2 * (sum(ptm_generator_mat[m]) - sum(ptm_generator_mat[m]; dims=2)[i])
                )
                for j in 1:3
                    ptm_parameter_mat[m][i, j] = exp(
                        -2 * (
                            sum(ptm_generator_mat[m]; dims=2)[i] +
                            sum(ptm_generator_mat[m]; dims=1)[j] -
                            2 * ptm_generator_mat[m][i, j]
                        ),
                    )
                end
            end
        end
        return ptm_parameter_vec, ptm_parameter_mat
    else
        error("Uneven lenghts in parameter arrays")
    end
end

function noise_ptm_coefficients(model::SPLNoiseModel)
    return noise_ptm_coefficients(noise_ptm_generators(model)...)
end

function noise1qblocks(sites::Vector{<:Index}, ptm_coefficients_vec)
    N = length(ptm_coefficients_vec)  # This is the number of qbits
    # Put the ptm_coefficients_vec[m] vector on the diagonal.
    q1g_ptm_matrix = [diagm(0 => ptm_coefficients_vec[m]) for m in 1:N]
    return [ITensor(q1g_ptm_matrix[m], sites[m]', dag(sites[m])) for m in 1:N]
end

function noise2qblocks(sites::Vector{<:Index}, ptm_coefficients_mat)
    N = length(ptm_coefficients_mat)  # This is the number of qbits, minus one
    # Unravel the ptm_coefficients_mat[m] matrix by rows, i.e. obtain the vector
    #
    #   ⎛ ptm_coefficients_mat[m][0,0] ⎞
    #   ⎜ ptm_coefficients_mat[m][0,1] ⎟
    #   ⎜ ptm_coefficients_mat[m][0,2] ⎟
    #   ⎜ ptm_coefficients_mat[m][0,3] ⎟
    #   ⎜ ptm_coefficients_mat[m][1,0] ⎟
    #   ⎝            ...               ⎠
    #
    # and put it on the diagonal.
    q2g_ptm_matrix = [diagm(0 => [ptm_coefficients_mat[m]'...]) for m in 1:N]
    return [
        ITensor(
            q2g_ptm_matrix[m], sites[m + 1]', sites[m]', dag(sites[m + 1]), dag(sites[m])
        ) for m in 1:N
    ]
end

"""
    noiselayer(sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat)

Return the MPO form of a layer of the noise map generated by the parameter sets
`ptm_generator_vec` and `ptm_generator_mat`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function noiselayer(sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat)
    if !(length(ptm_generator_vec) == length(ptm_generator_mat) + 1 == length(sites))
        error("Uneven lenghts in parameter arrays")
    end
    nqbits = length(ptm_generator_vec)

    ptm_coefficients_vec, ptm_coefficients_mat = noise_ptm_coefficients(
        ptm_generator_vec, ptm_generator_mat
    )

    q1g = noise1qblocks(sites, ptm_coefficients_vec)  # Single-qbit gates
    q2g = noise2qblocks(sites, ptm_coefficients_mat)  # Two-qbit gates

    noise = MPO(q1g)
    # This is just single-qbit gates. Now we apply all the odd gates, then
    # the even ones.
    for m in 1:2:(nqbits - 1)
        noise = apply(q2g[m], noise)
    end
    for m in 2:2:(nqbits - 1)
        noise = apply(q2g[m], noise)
    end

    # We know from the cited paper that the expected bond dimension of the MPO is 4,
    # so we happily truncate it without fear of any consequence (right??).
    return truncate(noise; maxdim=4)
end

"""
    noiselayer(sites::Vector{<:Index}, model::SPLNoiseModel)

Return the MPO form of a layer of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function noiselayer(sites::Vector{<:Index}, model::SPLNoiseModel)
    return noiselayer(sites, noise_ptm_generators(model)...)
end

# Naive MPO construction
# ----------------------
"""
    inversenoiselayer(sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat)

Return the MPO form of the inverse of the noise map generated by the parameter sets
`ptm_generator_vec` and `ptm_generator_mat`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function inversenoiselayer(sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat)
    if !(length(ptm_generator_vec) == length(ptm_generator_mat) + 1 == length(sites))
        error("Uneven lenghts in parameter arrays")
    end
    nqbits = length(ptm_generator_vec)

    ptm_coefficients_vec, ptm_coefficients_mat = noise_ptm_coefficients(
        -ptm_generator_vec, -ptm_generator_mat
    )
    # We obtain the inverse of the noise operators by supplying the same parameters of the
    # noise model but with the opposite sign.

    q1g = noise1qblocks(sites, ptm_coefficients_vec)  # Single-qbit gates
    q2g = noise2qblocks(sites, ptm_coefficients_mat)  # Two-qbit gates

    noiseinverse = MPO(q1g)
    # This is made just with single-qbit gates. Now we apply all the odd gates, then
    # the even ones.
    for m in 1:2:(nqbits - 1)
        noiseinverse = apply(q2g[m], noiseinverse)
    end
    for m in 2:2:(nqbits - 1)
        noiseinverse = apply(q2g[m], noiseinverse)
    end

    # We know from the cited paper that the expected bond dimension of the MPO is 4,
    # so we happily truncate it without fear of any consequence (right??).
    return truncate(noiseinverse; maxdim=4)
end

"""
    inversenoiselayer(sites::Vector{<:Index}, model::SPLNoiseModel)

Return the MPO form of the inverse of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function inversenoiselayer(sites::Vector{<:Index}, model::SPLNoiseModel)
    return inversenoiselayer(sites, noise_ptm_generators(model)...)
end

# Now we build the complete TEM MPO.
# We need the circuit (i.e. the gate list) and the noise MPO.
# The idea is that each gate comes with some noise.

"""
    tem_mpo(code::AbstractString,noiseptm_generator_vec,noiseptm_generator_mat;kwargs...)

Return the tensor-network error-mitigation MPO starting from the circuit structure
described in `code` (which contains the text of an OpenQASM source code file) and from
a sparse Pauli-Lindblad noise model generated by the parameters `noiseptm_generator_vec` and `noiseptm_generator_mat`.
Keyword arguments are forwarded to the `apply` function, so that it is possible to specify
the parameters for the contraction sequence (i.e. cutoff, maximum bond dimension).
"""
function tem_mpo(
    code::AbstractString, noiseptm_generator_vec, noiseptm_generator_mat; kwargs...
)
    sites, layers = gate_layers(code, "vQubit")
    # Each layer is made of some ITensors, which do not necessarily cover the whole width
    # of the circuit, so we call our `fullMPO` instead of the standard MPO constructor.
    layer_mpos = [fullMPO(sites, layer) for layer in layers]
    𝓝⁻¹ = inversenoiselayer(sites, noiseptm_generator_vec, noiseptm_generator_mat)

    # Middle-out contraction sequence
    # -------------------------------
    # M[k] = U[k] M[k-1] U[k]^-1 N[k]^-1,
    # M[0] = Id
    #
    # U[k] is the k-th layer of our circuit, so basically layer_mpos[k], while N[k]^-1
    # is our inverse noise MPO built following the sparse-Pauli-Lindblad model, for all k.
    tem = MPO(sites, "Id")
    for 𝓤 in layers
        tem = apply(𝓤, tem; apply_dag=true, kwargs...)
        tem = apply(tem, 𝓝⁻¹; kwargs...)
    end
    return tem
end

"""
    havecommonelements(A, B)

Return `true` if `A` and `B` have at least an element in common.
"""
function havecommonelements(A, B)
    for a in A
        if a in B
            return true
        end
    end
    return false
end

"""
    fullMPO(sites::Array{<:Index}, tensors::Vector{ITensor})

Create an MPO defined on the whole `sites` list, containing the given `tensors`.

The MPO is extended by filling the empty spots with an "Id" operator, therefore
an operator with OpName "Id" is required to be defined for the SiteTypes of the
remaining sites.

# Arguments
- `sites::Array{<:Index}`: the sites of the whole system.
- `tensors::Vector{ITensor}`: a list of tensors, with no common indices.
"""
function fullMPO(sites::Array{<:Index}, tensors::Vector{ITensor})
    # Compute the sites on which the MPO is defined
    tensors_sites = collect(Iterators.flatten(inds.(tensors; plev=0)))
    sites_not_in_tensors = findall(!in(tensors_sites), sites)
    # ↖ returns the indices (within the `sites` vector) of those sites which aren't
    #   in `tensor_sites`

    if !issubset(tensors_sites, sites)
        # Throw an error if some Index in one of the `tensors` is not in `sites`
        throw(BoundsError(tensors_sites, sites))
    end

    # Create an MPO filling the empty spots with an identity operator. The order of the
    # tensors should not matter, so we just put them at the end of the list.
    # Ideally we would just return
    #   MPO([tensors; [op("Id", sites, n) for n in sites_not_in_tensors]]
    # but if all sites are covered by `tensors` then the second list is empty; to be
    # precise, it is equal to a vector Any[]. The MPO constructor then rejects the argument
    # since the whole thing is a Vector{Any}. For this reason, we must treat this case
    # separately.
    if issetequal(tensors_sites, sites)  # They might not be in the same order!
        return MPO(tensors)
    else
        return MPO([tensors; [op("Id", sites, n) for n in sites_not_in_tensors]])
    end
end
