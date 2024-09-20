function pauli_string_constructors()
    # Test that the different PauliString constructors give the same result.
    len = 100
    nfactors = 10
    chars = rand(['X', 'Y', 'Z'], nfactors)
    positions = unique(sort(rand(1:100, nfactors)))

    # Compact form
    str1 = join(Iterators.flatten(zip(chars, string.(positions))))

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
