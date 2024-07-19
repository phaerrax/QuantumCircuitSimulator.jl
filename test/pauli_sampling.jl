function pauli_sampling(; nsamples=100)
    # Construct a random Pauli string and sample from the MPS built from it.
    # We expect to find always the same string in the samples, with overlap 1.
    nsites = 10
    sites = siteinds("vQubit", nsites)
    p = PauliString(mod.(rand(UInt8, nsites), 4))
    xp = MPS(sites, p)
    sampled_ps, overlaps = samplepaulistrings(xp, nsamples)
    return all(==(p), sampled_ps) && all(≈(1), overlaps)
end
