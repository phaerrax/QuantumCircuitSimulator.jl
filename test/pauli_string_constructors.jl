function pauli_string_constructors()
    # Test that the different PauliString constructors give the same result.
    len = 100
    chars = unique(sort(rand(['I', 'X', 'Y', 'Z'], 10)))
    nfactors = length(chars)
    positions = rand(1:100, nfactors)

    # Compact form
    str1 = reduce(*, collect(Iterators.flatten(zip(chars, string.(positions)))))

    # Full form
    str2_factors = repeat(['I'], len)
    for (pos, char) in zip(positions, chars)
        str2_factors[pos] = char
    end
    str2 = *(str2_factors...)

    return PauliString(len, zip(chars, positions)...) ==
           PauliString(len, str1) ==
           PauliString(str2)
end
