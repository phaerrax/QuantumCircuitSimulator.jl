"""
    SPLNoiseModel

A sparse Pauli-Lindblad noise model, represented by a dictionary of generators and a list
of qbits it acts on.
"""
struct SPLNoiseModel
    parameters::Dict{PauliString,Real}
    siteindices::Vector{Int}
    # For each Pauli string `p` in `keys(parameters)`, the factor `p[j]` acts on the qbit
    # at position `siteindices[j]`.
    function SPLNoiseModel(dict::Dict{PauliString,<:Real}, indices)
        if !allequal(length.(keys(dict)))
            # Check that all strings have the same length...
            error("Pauli strings don't have the same length")
        elseif !all(
            n -> (n == 1 || n == 2), [PauliStringTensors.order(p.first) for p in dict]
        )
            # ...and that their order is 1 or 2
            error("Some Pauli strings don't have order 1 or 2")
        elseif length(first(keys(dict))) != length(indices)
            error("Length of Pauli strings and number of indices don't match")
        else
            new(dict, indices)
        end
    end
end

function SPLNoiseModel(dict::Dict{PauliString,<:Real})
    strlen = length(first(keys(dict)))
    return SPLNoiseModel(dict, 1:strlen)
end

Base.length(model::SPLNoiseModel) = length(model.siteindices)
nqbits(model::SPLNoiseModel) = length(model.siteindices)
qbitsites(model::SPLNoiseModel) = model.siteindices

"""
    crop(model::SPLNoiseModel, range)

Return a subset of `model` keeping only the Pauli strings whose non-trivial factors lie
completely within `range`.
"""
function crop(model::SPLNoiseModel, range)
    if !issubset(range, qbitsites(model))
        error("New range not contained in the domain of the noise model.")
    end
    cropped = Dict{PauliString,Real}()
    for (p, v) in model.parameters
        nontrivialinds = indices(p)
        if issubset(nontrivialinds, range)
            push!(cropped, crop(p, range) => v)
        end
    end
    return SPLNoiseModel(cropped, range)
end

function SPLNoiseModel(file::AbstractString)
    dict = JSON.parsefile(file)
    return SPLNoiseModel(Dict(PauliString(str) => coeff for (str, coeff) in dict))
end

function SPLNoiseModel(file::AbstractString, indices)
    dict = JSON.parsefile(file)
    return SPLNoiseModel(Dict(PauliString(str) => coeff for (str, coeff) in dict), indices)
end

function noise_ptm_generators(model::SPLNoiseModel)
    N = nqbits(model)

    vec = [zeros(Float64, 3) for _ in 1:N]
    mat = [zeros(Float64, 3, 3) for _ in 1:(N - 1)]
    for (k, v) in model.parameters
        if PauliStringTensors.order(k) == 1
            # indices(k): index of non-trivial factors in Pauli string
            # operators(k): non-trivial factors in Pauli string
            vec[first(indices(k))][operators(k)...] = v
        elseif PauliStringTensors.order(k) == 2
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
    N = length(ptm_coefficients_vec)  # number of qbits
    # Put the ptm_coefficients_vec[m] vector on the diagonal.
    q1g_ptm_matrix = [diagm(0 => ptm_coefficients_vec[m]) for m in 1:N]
    return [ITensor(q1g_ptm_matrix[m], sites[m]', dag(sites[m])) for m in 1:N]
end

function noise2qblocks(sites::Vector{<:Index}, ptm_coefficients_mat)
    N = length(ptm_coefficients_mat)  # number of qbits, minus one
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
    noiselayer(
        sites::Vector{<:Index},
        indices,
        ptm_generator_vec,
        ptm_generator_mat;
        maxbonddim=nothing
    )

Return the MPO form of a layer of the noise map generated by the parameter sets
`ptm_generator_vec` and `ptm_generator_mat`, as explained in [1], acting on the indices
`sites[i]` for all `i` in `indices`. The MPO will be defined on all `sites` (as the
identity on those sites `i` which are not in `indices`).
Optionally, truncate the resulting MPO to a bond dimension equal to `maxbonddim`.

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function noiselayer(
    sites::Vector{<:Index},
    indices,
    ptm_generator_vec,
    ptm_generator_mat;
    maxbonddim=nothing,
)
    if length(ptm_generator_vec) != length(ptm_generator_mat) + 1 ||
        length(ptm_generator_vec) > length(indices) ||
        length(ptm_generator_mat) + 1 > length(indices)
        error("Incompatible lenghts in parameter arrays")
    end
    nqbits = length(ptm_generator_vec)

    ptm_coefficients_vec, ptm_coefficients_mat = noise_ptm_coefficients(
        ptm_generator_vec, ptm_generator_mat
    )

    q1g = noise1qblocks(sites[indices], ptm_coefficients_vec)  # Single-qbit gates
    q2g = noise2qblocks(sites[indices], ptm_coefficients_mat)  # Two-qbit gates

    noise = MPO(sites, "Id")
    for m in 1:nqbits
        noise = apply(q1g[m], noise)
    end
    # Now we apply all the odd gates, then the even ones.
    for m in 1:2:(nqbits - 1)
        noise = apply(q2g[m], noise)
    end
    for m in 2:2:(nqbits - 1)
        noise = apply(q2g[m], noise)
    end

    if !isnothing(maxbonddim) && maxbonddim > 0
        truncate!(noise; maxdim=maxbonddim)
    end

    return noise
end

"""
    noiselayer(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)

Return the MPO form of a layer of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`, as explained in [1].
Optionally, truncate the resulting MPO to a bond dimension equal to `maxbonddim`.

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function noiselayer(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return noiselayer(
        sites, qbitsites(model), noise_ptm_generators(model)...; maxbonddim=maxbonddim
    )
end

"""
    MPO(sites::Vector{<:Index}, noisemodel::SPLNoiseModel; maxbonddim)

Return the MPO form of a layer of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`, as explained in [1].
Optionally, truncate the resulting MPO to a bond dimension equal to `maxbonddim`.

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function ITensorMPS.MPO(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return noiselayer(sites, model; maxbonddim=maxbonddim)
end

# Naive MPO construction
# ----------------------
"""
    inversenoiselayer(
        sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat; maxbonddim=nothing
    )

Return the MPO form of the inverse of the noise map generated by the parameter sets
`ptm_generator_vec` and `ptm_generator_mat`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function inversenoiselayer(
    sites::Vector{<:Index},
    indices,
    ptm_generator_vec,
    ptm_generator_mat;
    maxbonddim=nothing,
)
    # TODO we can reuse `noiselayer` here, if we just feed it the opposites of
    # the generator vectors.
    if length(ptm_generator_vec) != length(ptm_generator_mat) + 1 ||
        length(ptm_generator_vec) > length(indices) ||
        length(ptm_generator_mat) + 1 > length(indices)
        error("Incompatible lenghts in parameter arrays")
    end
    nqbits = length(ptm_generator_vec)

    ptm_coefficients_vec, ptm_coefficients_mat = noise_ptm_coefficients(
        -ptm_generator_vec, -ptm_generator_mat
    )
    # We obtain the inverse of the noise operators by supplying the same parameters of the
    # noise model but with the opposite sign.

    q1g = noise1qblocks(sites[indices], ptm_coefficients_vec)  # Single-qbit gates
    q2g = noise2qblocks(sites[indices], ptm_coefficients_mat)  # Two-qbit gates

    noiseinverse = MPO(sites, "Id")
    for m in 1:nqbits
        noiseinverse = apply(q1g[m], noiseinverse)
    end
    # Now we apply all the odd gates, then the even ones.
    for m in 1:2:(nqbits - 1)
        noiseinverse = apply(q2g[m], noiseinverse)
    end
    for m in 2:2:(nqbits - 1)
        noiseinverse = apply(q2g[m], noiseinverse)
    end

    if !isnothing(maxbonddim) && maxbonddim > 0
        truncate!(noiseinverse; maxdim=maxbonddim)
    end

    return noiseinverse
end

"""
    inversenoiselayer(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)

Return the MPO form of the inverse of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function inversenoiselayer(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return inversenoiselayer(
        sites, qbitsites(model), noise_ptm_generators(model)...; maxbonddim=maxbonddim
    )
end

"""
    inverseMPO(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)

Return the MPO form of the inverse of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`, as explained in [1].

# References

[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function inverseMPO(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return inversenoiselayer(sites, model; maxbonddim=maxbonddim)
end

# Now we build the complete TEM MPO.
# We need the circuit (i.e. the gate list) and the noise MPO.
# The idea is that each gate comes with some noise.

"""
    tem_mpo(
    code::OpenQASM.Types.MainProgram,
    noiseptm_generator_vec,
    noiseptm_generator_mat;
    kwargs...,
)

Return the tensor-network error-mitigation MPO starting from the circuit structure
described in `code` (which contains the text of an OpenQASM source code file) and from
a sparse Pauli-Lindblad noise model generated by the parameters `noiseptm_generator_vec` and `noiseptm_generator_mat`.
Keyword arguments are forwarded to the `apply` function, so that it is possible to specify
the parameters for the contraction sequence (i.e. cutoff, maximum bond dimension).
"""
function tem_mpo(
    code::OpenQASM.Types.MainProgram,
    noiseptm_generator_vec,
    noiseptm_generator_mat;
    kwargs...,
)
    sites, unitarylayers = gatelayers(code, "vQubit")

    # Build the full noise layer with the appropriate function.
    inversenoise = inversenoiselayer(sites, noiseptm_generator_vec, noiseptm_generator_mat)
    # TODO add truncation parameters here ↗

    # Middle-out contraction sequence
    # -------------------------------
    # M[k] = U[k] M[k-1] U[k]^-1 N[k]^-1,
    # M[0] = Id
    #
    # U[k] is the k-th layer of our circuit, which is stored in 𝓤[k] as a vector of ITensor
    # objects. We don't need to create an MPO out of them, we can just apply them directly
    # one by one to the "main" MPO (they all commute so we don't need to care about the
    # order of composition); N[k]^-1 is our inverse-noise MPO built following the
    # sparse-Pauli-Lindblad model (the same for all k).
    tem = MPO(sites, "Id")
    for ul in unitarylayers
        tem = apply(ul, tem; apply_dag=true, kwargs...)
        tem = apply(tem, inversenoise; kwargs...)
    end
    return tem
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
