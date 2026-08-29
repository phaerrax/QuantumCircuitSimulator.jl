"""
    SPLNoiseModel

A sparse Pauli-Lindblad noise model, represented by a dictionary of generators and a list
of qbits it acts on.
"""
struct SPLNoiseModel
    parameters::Dict{PauliString,<:Real}
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

function SPLNoiseModel(
    dict::Dict{<:AbstractString,<:Real}, indices=1:length(first(keys(dict)))
)
    dict = Dict(PauliString(p) => coeff for (p, coeff) in dict)
    return SPLNoiseModel(dict, indices)
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
`ptm_generator_vec` and `ptm_generator_mat`, acting on the indices `sites[i]` for all `i` in
`indices`. The MPO will be defined on all `sites` (as the identity on those sites `i` which
are not in `indices`).
Optionally, truncate the resulting MPO to a bond dimension equal to `maxbonddim`.

# References

* [Filippov2023](@cite) S. Filippov, M. Leahy, M. A. C. Rossi and G. García-Pérez, _arXiv_
  2307.11740 (2023)
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
described in the sparse Pauli-Lindblad model object `model`.
Optionally truncate the resulting MPO to a bond dimension equal to `maxbonddim`.
"""
function noiselayer(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return noiselayer(
        sites, qbitsites(model), noise_ptm_generators(model)...; maxbonddim=maxbonddim
    )
end

"""
    MPO(sites::Vector{<:Index}, noisemodel::SPLNoiseModel; maxbonddim)

Return the MPO form of a layer of the noise map generated by the parameter sets
described in the sparse Pauli-Lindblad model object `model`.
Optionally truncate the resulting MPO to a bond dimension equal to `maxbonddim`.
"""
function ITensorMPS.MPO(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return noiselayer(sites, model; maxbonddim=maxbonddim)
end

# Inverse noise layer construction

"""
    inversenoiselayer(
        sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat; maxbonddim=nothing
    )

Return the MPO form of the inverse of the noise map generated by the parameter sets
`ptm_generator_vec` and `ptm_generator_mat`.

# References

* [Filippov2023](@cite) S. Filippov, M. Leahy, M. A. C. Rossi and G. García-Pérez, _arXiv_
  2307.11740 (2023)
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
described in the sparse Pauli-Lindblad model object `model`.
"""
function inversenoiselayer(sites::Vector{<:Index}, model::SPLNoiseModel; maxbonddim=nothing)
    return inversenoiselayer(
        sites, qbitsites(model), noise_ptm_generators(model)...; maxbonddim=maxbonddim
    )
end
