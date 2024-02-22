using ITensors, PseudomodesTTEDOPA

"""
    ITensors.space(st::SiteType"vQubit"; dim = 2)

Create the Hilbert space for a site of type "vQubit", i.e. a vectorised
qbit, where the vectorisation is performed wrt the generalised
Gell-Mann basis of `Mat(ℂ²)`, composed of Hermitian traceless matrices
together with the identity matrix.
"""
ITensors.space(::SiteType"vQubit") = 4

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vQubit")
    v = ITensors.state(sn, SiteType("Qubit"))
    return PseudomodesTTEDOPA.vec(kron(v, v'), gellmannbasis(2))
end
function vop(sn::StateName, ::SiteType"vQubit")
    sn = statenamestring(sn)
    on = sn[1] == 'v' ? sn[2:end] : sn
    return PseudomodesTTEDOPA.vec(PseudomodesTTEDOPA.try_op(OpName(on), SiteType("Qubit")), gellmannbasis(2))
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"0", st::SiteType"vQubit") = vstate(sn, st)
ITensors.state(sn::StateName"1", st::SiteType"vQubit") = vstate(sn, st)

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vQubit")
    return PseudomodesTTEDOPA.vec(x -> mat * x, gellmannbasis(2))
end
function postmultiply(mat, ::SiteType"vQubit")
    return PseudomodesTTEDOPA.vec(x -> x * mat, gellmannbasis(2))
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the S=1/2 site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vQubit"; kwargs...)
    name = strip(String(ITensors.name(on))) # Remove extra whitespace
    if name == "Id"
        return Matrix(1.0I, 4, 4)
    end
    dotloc = findfirst("⋅", name)
    # This returns the position of the cdot in the operator name String.
    # It is `nothing` if no cdot is found.
    if !isnothing(dotloc)
        on1, on2 = nothing, nothing
        on1 = name[1:prevind(name, dotloc.start)]
        on2 = name[nextind(name, dotloc.start):end]
        # If the OpName `on` is written correctly, i.e. it is either "A⋅" or "⋅A" for some
        # A, then either `on1` or `on2` has to be empty (not both, not neither of them).
        if (on1 == "" && on2 == "") || (on1 != "" && on2 != "")
            throw(
                ArgumentError(
                    "Invalid operator name: $name. Operator name is not \"Id\" or of the " *
                    "form \"A⋅\" or \"⋅A\"",
                ),
            )
        end
        # name == "⋅A" -> on1 is an empty string
        # name == "A⋅" -> on2 is an empty string
        if on1 == ""
            mat = PseudomodesTTEDOPA.try_op(OpName(on2), SiteType("S=1/2"); kwargs...)
            return postmultiply(mat, st)
        elseif on2 == ""
            mat = PseudomodesTTEDOPA.try_op(OpName(on1), SiteType("S=1/2"); kwargs...)
            return premultiply(mat, st)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
end

#=
Some gates require parameters. We need to pass these parameters to gkslcommutator_itensor
as well, if we want to build the commutator.

Example 1: a commutator with the U gate
---------------------------------------
We need to get to the expression (with some random values for the angles)
    -im * (op("U⋅", s[n]; θ=pi/4, ϕ=0, λ=pi/8) - op("⋅U", s[n]; θ=pi/4, ϕ=0, λ=pi/8)).
The parameters are best associated to the gate name, rather than to the site, therefore
it makes sense to redefine gkslcommutator_itensor so that it accepts arguments as
    gkslcommutator_itensor(sites, X, n, X', n'...)
where X, X',... contain information about the gate name _and_ its parameters.
Since the parameters are given to ITensors.op as keyword arguments, i.e. a NamedTuple,
it might make sense to specify X as a Tuple with the operator name (a String) in the first
position and the parameters (a NamedTuple) in the second position:
    X = ("U", (θ=pi/4, ϕ=0, λ=pi/8))
and later build the relevant operators using
    on, args = X
    lmult = ITensors.op("$on⋅", s, n; args...)
    rmult = ITensors.op("⋅$on", s, n; args...)
    -im * (lmult - rmult)

Example 2: a commutator with the H gate
---------------------------------------
What if the gate doesn't have any parameter? Well, luckily the above syntax works even if
`args` is empty, but we still need to provide it, so we will use
    X = ("H", ())
and then go on as above; `args` will be an empty tuple, and we have
    ITensors.op("H⋅", s, n'; ()...) == ITensors.op("H⋅", s, n').

Example 3: both the above gates at the same time
------------------------------------------------
Note that the `gkslcommutator_itensor(::Vector{<:Index}, ::Tuple{String,Int}...)` function
(the original one) works like this: you pass the vector of site indices and then some
of {String,Int} pairs, which get slurped into a Vector of Pair{String,Int} items:
    gkslcommutator_itensor(s, ("A", 1), ("B", 3)).
For the purposes of this library, "A" and "B" need to be replaced by the X and X' things
above, so by (another!) tuple:
    gkslcommutator_itensor(s, (("U", (θ=pi/4, ϕ=0, λ=pi/8)), 1), (("H", ()), 3)).
The second and third argument of `gkslcommutator_itensor` here (try this on the REPL!) are
Tuple{Tuple{String, Any}, Int} items. The new signature should maybe then be
    gkslcommutator_itensor(::Vector{<:Index}, ::Tuple{Tuple{String, Any},Int}...).
=#

function gkslcommutator_itensor(sites::Vector{<:Index}, items::Tuple{Tuple{String, Any},Int}...)
    # Unpacking the arguments:
    # - separate operator names/args and site numbers
    operators = first.(items)
    site_numbers = last.(items)
    if !allunique(site_numbers)
        error("Some sites are repeated in the list. Please use unique site indices.")
        # This possibility is not allowed for now, since it would create issues in the
        # multiplication loop below. Basically, if two ITensors with the same indices
        # are multiplied together (with a simple `*`) then both indices get contracted
        # and we end up with a scalar. We should use `apply` instead, but it does not
        # work with OneITensor objects...
    end
    # - separate operator names from parameters
    operator_names = first.(operators)
    operator_args = last.(operators)

    lmult = ITensors.OneITensor()
    rmult = ITensors.OneITensor()
    for (on, args, j) in zip(operator_names, operator_args, site_numbers)
        lmult *= op("$on⋅", sites, j; args...)
        rmult *= op("⋅$on", sites, j; args...)
    end
    return -im * (lmult - rmult)
end

#function gkslcommutator_itensor(sites::Vector{<:Index}, args...)
#    return gkslcommutator_itensor(sites, makeopsumpairs(args...)...)
#end
