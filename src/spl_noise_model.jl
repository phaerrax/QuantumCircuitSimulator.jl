"""
    noise_ptm_coefficients(ptm_generator_vec, ptm_generator_mat)

Return the parameters of the noise model in the sparse Pauli-Lindblad model from the
given set of coefficients `ptm_generator_vec` and `ptm_generator_mat`, respectively a set
of ``N`` 3-element vectors and a set of ``N`` 3x3 matrices, where ``N`` represents the
number of qbits in the circuit.
"""
function noise_ptm_coefficients(ptm_generator_vec, ptm_generator_mat)
    if length(ptm_generator_vec) == size(ptm_generator_mat, 1) == size(ptm_generator_mat, 2)
        n_qbits = length(ptm_generator_vec)
        ptm_parameter_vec = [OffsetArray(ones(4), 0:3) for _ in 1:n_qbits]
        ptm_parameter_mat = [OffsetArray(ones(4, 4), 0:3, 0:3) for _ in 1:n_qbits]
        for m in 1:n_qbits
            for i in 1:3
                ptm_parameter_vec[m][i] = exp(
                    -2 * sum([ptm_generator_vec[m][j] for j in 1:3 if j != i])
                )
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

# Naive MPO construction
# ----------------------
"""
    noise_inverse_mpo(sites::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat)

Return the MPO form of the inverse of the noise map generated by the parameter sets
`ptm_generator_vec` and `ptm_generator_mat`, as explained in [1].

# References
[1] 'Scalable tensor-network error mitigation for near-term quantum computing',
Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, and Guillermo Garcı́a-Pérez.
"""
function noise_inverse_mpo(s::Vector{<:Index}, ptm_generator_vec, ptm_generator_mat)
    if !(length(ptm_generator_vec) == length(ptm_generator_mat) == length(s))
        error("Uneven lenghts in parameter arrays")
    end
    n_qbits = length(ptm_generator_vec)

    ptm_parameter_vec, ptm_parameter_mat = noise_ptm_coefficients(
        ptm_generator_vec, ptm_generator_mat
    )

    # Single-qbit gates
    # -----------------
    # Put the ptm_parameter_vec[m] vector on the diagonal.
    q1g_ptm_matrix = [diagm(0 => ptm_parameter_vec[m]) for m in 1:n_qbits]
    q1g = [ITensor(q1g_ptm_matrix[m], s[m]', dag(s[m])) for m in 1:n_qbits]

    # Two-qbit gates
    # --------------
    # Unravel the ptm_parameter_mat[m] matrix by rows, i.e. obtain the vector
    #
    #   ⎛ ptm_parameter_mat[m][0,0] ⎞
    #   ⎜ ptm_parameter_mat[m][0,1] ⎟
    #   ⎜ ptm_parameter_mat[m][0,2] ⎟
    #   ⎜ ptm_parameter_mat[m][0,3] ⎟
    #   ⎜ ptm_parameter_mat[m][1,0] ⎟
    #   ⎝         ...               ⎠
    #
    # and put it on the diagonal.
    q2g_ptm_matrix = [diagm(0 => [ptm_parameter_mat[m]'...]) for m in 1:n_qbits]
    q2g = [
        ITensor(q2g_ptm_matrix[m], s[m + 1]', s[m]', dag(s[m + 1]), dag(s[m])) for
        m in 1:(n_qbits - 1)
    ]

    invN = MPO(q1g)
    # This is made just with single-qbit gates. Now we apply all the odd gates, then
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

# Now we build the complete TEM MPO.
# We need the circuit (i.e. the gate list) and the noise MPO.
# The idea is that each gate comes with some noise, but how exactly do we layer the unitary
# gates and the noise MPO? We could (but I don't know if this is right) group in a single
# layer all subsequent gates that commute with each other; they can be put "in a row", so
# to speak.
# We can detect this by looking at the Indices of the gate ITensors: as soon as we see a
# repeated index, we stop grouping them and insert a noise layer.
# (Circuit barriers could come into play at this point...)

# First, group the gates into layers.
function gate_layers(code::AbstractString, st::AbstractString)
    sites, gate_list = gates(code, st)

    # We put the gates in a stack, so that we can pop them one at a time. We reverse the
    # order, so that the first gate in the list is also the first in the stack.
    gate_stack = Stack{ITensor}()
    foreach(g -> push!(gate_stack, g), reverse(gate_list))

    layers = Vector{ITensor}[]
    # Will contain only gates, no noise MPOs. Each layer will be a vector of gates.

    # Create a new empty layer.
    current_layer = ITensor[]
    while !isempty(gate_stack)
        # Look up the next tensor in the gate sequence: it is `first(gate_stack)`.
        # Does this gate have an Index which is already present in the current layer?
        # 1. inds(first(gate_stack)) is a tuple of Index objects
        # 2. inds.(current_layer) is a _vector_ of tuples of Index objects
        # we need both of them to be a big list of Index, not of tuples:
        # 1. `collect` transforms a Tuple into a Vector
        # 2. `flatten` transforms a Vector of Tuple{X,X,...} into a Vector{X}
        if havecommonelements(
            collect(inds(first(gate_stack))), Iterators.flatten(inds.(current_layer))
        )
            # If yes: stop adding the layer, save it and start a new layer with this gate
            push!(layers, current_layer)
            current_layer = [pop!(gate_stack)]
        else
            # If not: add the gate to the current layer and go on.
            push!(current_layer, pop!(gate_stack))
        end
    end

    return sites, layers
end

# Then, build the MPO interleaving gates and noise layers.
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
    𝓝⁻¹ = noise_inverse_mpo(sites, noiseptm_generator_vec, noiseptm_generator_mat)

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