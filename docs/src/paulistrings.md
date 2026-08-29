# Pauli-string sampling

When an MPS represents a vectorised density matrix or and observable in the Pauli
transfer matrix (PTM) basis, it can be useful to inspect it in terms of its
Pauli-string components, rather than as a tensor network.
The [`samplepaulistrings`](@ref) function draws Pauli strings from such an MPS
(interpreting it as a probability distribution over strings) and returns, for
each sampled string, its coefficient in the expansion of the MPS as a linear
combination of Pauli strings:

```jldoctest paulistrings
julia> using ITensors, ITensorMPS, QuantumCircuitSimulator

julia> sites = siteinds("vQubit", 4);

julia> v = random_mps(sites);

julia> strings, coefficients = samplepaulistrings(v, 1000);
```

The resulting `coefficients[k]` is the coefficient of `strings[k]` in the
decomposition of `v`.

Building on top of this result, the [`relevantpaulistrings`](@ref) function
returns a more directly useful summary: it samples a set of Pauli strings,
computes the frequency with which they appear, and returns a list of `(string,
coefficient, frequency)` tuples sorted from the largest to the smallest
coefficient (in absolute value).  It can optionally be capped to a maximum
number of strings (`maxn`) or filtered by a minimum sampling frequency
(`cutoff`):

```jldoctest paulistrings
julia> components = relevantpaulistrings(v; nsamples=1000, maxn=10);
```

This is typically used to obtain a compact, human-readable approximation of an
otherwise intractably large observable or density matrix, keeping only its most
significant Pauli-string components.
