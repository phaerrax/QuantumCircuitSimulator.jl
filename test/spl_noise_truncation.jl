function random_two_site_paulistring(N)
    ps = zeros(UInt8, N)
    lind = rand(1:N)
    ps[lind] = rand(1:3)  # first site is always ≠ I
    if lind < N
        ps[lind + 1] = rand(0:3)
        # second site can be I: this way we get both 1-site and 2-site Pauli strings.
    end
    return PauliString(ps)
end

function test_spl_noise_truncation(; N=8, atol=1e-10)
    # We build an artificial SPL noise model with random strings and coefficients.
    spl_dictionary = Dict{PauliString,Float64}()
    # There are at most 3N + 9(N-1) distinct 1-site and (contiguous) 2-site Pauli strings
    # for N qbits. Let's draw half of this number for safety.
    max_n = 2N + 9(N - 1)
    for _ in 1:max_n
        push!(spl_dictionary, random_two_site_paulistring(N) => rand())
    end

    noisemodel = SPLNoiseModel(spl_dictionary)
    sites = siteinds("vQubit", nqbits(noisemodel))

    # We know that an SPL noise layer can be faithfully represented as an MPO with bond
    # dimension 4, but we check it anyway.
    noise_mpo = noiselayer(sites, noisemodel)
    truncated_noise_mpo = truncate(noise_mpo; maxdim=4)

    # The inverse noise MPO too can be exactly represented with bond dimension 4.
    invnoise_mpo = inversenoiselayer(sites, noisemodel)
    truncated_invnoise_mpo = truncate(invnoise_mpo; maxdim=4)

    identity = MPO(sites, "Id")

    @test noise_mpo ≈ truncated_noise_mpo
    @test invnoise_mpo ≈ truncated_invnoise_mpo
    @test apply(noise_mpo, invnoise_mpo) ≈ identity
    @test apply(invnoise_mpo, noise_mpo) ≈ identity

    # These two tests involve much more numerical work, so we need a looser tolerance.
    @test apply(truncated_noise_mpo, truncated_invnoise_mpo) ≈ identity atol=1e-8
    @test apply(truncated_invnoise_mpo, truncated_noise_mpo) ≈ identity atol=1e-8

    return nothing
end
