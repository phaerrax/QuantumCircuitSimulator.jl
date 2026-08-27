function test_noise_layer_preserves_trace(; N=50)
    v = Vector{Vector{Float64}}(undef, N)
    for i in 1:N
        v[i] = rand(3)
    end

    m = Vector{Matrix{Float64}}(undef, N - 1)
    for i in 1:(N - 1)
        m[i] = rand(3, 3)
    end

    sites = siteinds("vQubit", N)
    vid = MPS(sites, "Id")

    𝒩 = noiselayer(sites, 1:length(sites), v, m)

    x = random_mps(sites; linkdims=2)

    @test dot(vid, x) ≈ dot(vid, apply(𝒩, x))
end
