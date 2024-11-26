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

function _contractPTM(p::PauliString, v::MPS)
    # Note that this function does not return the inner product of the Pauli string and v.
    # In fact it computes the inner product of the e_{p_1,..., p_N} element of the PTM
    # basis with v, where N = length(v), and e_{p_1,..., p_N} == 2^(-N/2) * p.
    c = ITensors.OneITensor()
    for n in 1:length(v)
        c *= _contract(onehot(siteind(v, n) => p.string[n] + 1), v[n])
    end
    return scalar(c)
end

"""
    samplepaulistrings(v::MPS, nsamples::Integer)

Sample `nsamples` Pauli strings from `v` and compute their overlap with the MPS.
Return a pair `ps, overlaps` where `overlaps[k]` is the coefficient of the `ps[k]`
component of `v`: this means that if we write `v` as a linear combination of Pauli
strings ``v = ∑ₖ cₖσₖ`` then `overlaps[k]` is the coefficient ``cₖ``.

Note that the MPS of a Pauli string is not normalized in the Hilbert-Schmidt inner
product ``⟨A,B⟩ = tr(A† B)``: the norm of a Pauli string of length `N` is ``2^(N/2)``.
"""
function samplepaulistrings(v::MPS, nsamples::Integer)
    if any(t -> !(SiteType("vQubit") in t), _sitetypes.(siteinds(v)))
        error("samplepaulistrings works for vQubit site types only.")
    end

    # The sampling algorithm requires a normalized MPS orthogonalized on the first site.
    vn = orthogonalize(v, 1)
    vn /= norm(vn)
    sites = siteinds(v)

    ps = Vector{PauliString}(undef, nsamples)
    for i in 1:nsamples
        ps[i] = PauliString(sample(vn) .- 1)
        # sample(vn) gives us a vector of elements from {1, 2, 3, 4}; they are indices
        # referring to the PTM basis {I/sqrt(2), X/sqrt(2), Y/sqrt(2), Z/sqrt(2)}.
        # PauliString objects work with {0, 1, 2, 3} (respectively) instead so we need to
        # decrease by one.
    end

    # We compute the overlap of each Pauli string with the (normalized) MPS v.
    # Note that the MPS created from a PauliString object is not normalized, since it uses
    # the {I, X, Y, Z} matrices instead of their rescaled counterparts in the PTM basis.
    # For this reason, each MPS(sites, p) has norm 2^(length(sites)/2), and the overlap
    # must be normalized accordingly.
    # The orthonormal PTM basis is (N := length(sites))
    #   eₖ = 2^(-N/2) σₖ,
    # so we can extract the coefficients with the inner product as follows:
    #   v = ∑ₖ cₖσₖ
    #   ⟨eⱼ, v⟩ = ∑ₖ cₖ⟨eⱼ,σₖ⟩ = 2^(N/2) ∑ₖ cₖ⟨eⱼ,eₖ⟩ =  2^(N/2) cⱼ
    # therefore cₖ = 2^(-N/2) ⟨eₖ, v⟩ = 2^(-N) ⟨σₖ, v⟩.
    normf = 2^(length(sites) / 2)
    overlaps = Vector{ComplexF64}(undef, nsamples)
    for i in 1:nsamples
        overlaps[i] = _contractPTM(ps[i], v) / normf
    end

    return ps, overlaps
end

"""
    relevantpaulistrings(v::MPS; nsamples, maxn, cutoff, imag_atol=1e-10)

Sample at most `nsamples` Pauli strings from the MPS `v` (which is assumed to represent an
observable in the PTM basis) and return a list of tuples of the form `(ps, coeff, freq)`
where:

* `ps` is a Pauli string (a PauliString object)
* `coeff` is its coefficient within `v`
* `freq` is the frequency with which it was sampled (it should equal ``|coeff|²`` in the
    ``nsamples → ∞`` limit).

The list is shown from the most to least relevant component, i.e. highest to the lowest
modulus of the coefficient. It can be cutoff after a certain maximum number of
strings, or below a set frequency.
"""
function relevantpaulistrings(v::MPS; nsamples, maxn=0, cutoff=0, imag_atol=1e-10)
    # 1) Sample some Pauli strings from the observable
    strings, coefficients = samplepaulistrings(v, nsamples)
    # We can assume that the coefficients are real, and just print a warning if some of
    # them aren't.
    test_imag = findfirst(>(imag_atol), imag.(coefficients))
    if !isnothing(test_imag)
        @warn "Imaginary value above $imag_atol threshold: $test_imag"
    end
    ps = collect(zip(strings, real.(coefficients)))

    # 2) Count the occurrences of each (unique) string, and compute the frequencies
    # The function `count(t -> t[1] == p[1], ps)` returns how many elements in `ps` are such
    # that `t -> t[1] == p[1]` is true, i.e. how many of them have the string (which is the
    # first element in the tuple) equal to `p[1]`.
    # While we're at it, we convert the Pauli string to a String format.
    # (The sampling frequencies aren't actually _that_ useful, since we can efficiently
    # compute the coefficients anyway, but we keep them as a way to check the results.)
    frequencies = [
        (p[1], p[2], count(t -> t[1] == p[1], ps) / nsamples) for p in unique(first, ps)
    ]

    # 4) Order the strings by the absolute value of the coefficient (in decreasing order)
    sort!(frequencies; by=t -> abs(t[2]), rev=true)

    # 5) Keep only the first `maxn` strings, if specified
    mostfrequent = if maxn > 0
        first(frequencies, maxn)
    else
        frequencies
    end

    # 6) Remove strings that did not appear frequently enough
    cutoff > 0 && filter!(t -> last(t) > cutoff, mostfrequent)

    return mostfrequent
end

function samplepaulistrings_progress()
    @error "Please load the ProgressMeter package to load this function."
end

function relevantpaulistrings_progress()
    @error "Please load the ProgressMeter package to load this function."
end
