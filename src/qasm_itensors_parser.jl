using ITensors

"""
    splitlines(code::AbstractString)

Split code at semicolons, producing a list whose elements are single lines of code.

# Example

```julia-repl
julia> code = "OPENQASM 2.0;\nqreg a[2];\ncx a[0], a[1];\nz a[1];";

julia> splitlines(code)
4-element Vector{SubString{String}}:
 "OPENQASM 2.0"
 "qreg a[2]"
 "cx a[0], a[1]"
 "z a[1]"
```
"""
function splitlines(code::AbstractString)
    lines = split(code, ";")  # Split at semicolons
    lines = strip.(lines)  # Remove trailing and leading whitespace (such as stray '\n's)
    return filter(!isempty, lines)  # Remove empty lines
end

"""
    stripcomments(str::AbstractString)

Remove comment lines from the given code.

# Example
```julia-repl
julia> code = "OPENQASM 2.0;//this is a comment\nqreg a[2];\ncx a[0], a[1];\nz a[1];";

julia> TEM.stripcomments(code)
"OPENQASM 2.0;\nqreg a[2];\ncx a[0], a[1];\nz a[1];"
```
"""
function stripcomments(str::AbstractString)
    # Remove comments, i.e. everything between a '//' and a newline
    # (Strings are immutable, so we can't use `replace!` here)
    return replace(str, r"//[^\n]*" => "")
end

"""
    removepreprocessor!(lines::Vector{<:AbstractString})

Filter out from the list lines which are preprocessor instructions or file inclusions.
"""
function removepreprocessor!(lines::Vector{<:AbstractString})
    # Remove preprocessor lines
    filter!(s -> !occursin("OPENQASM 2.0", s), lines)
    return filter!(s -> !occursin("include", s), lines)
end

"""
    qbit_registers(code::AbstractString)

Return names and lengths of the quantum registers declared in the given
code.

# Example

```julia-repl
julia> code = "qreg a[3];\nqreg b[2];";

julia> qbit_registers(code)
(["a", "b"], [3, 2])
```
"""
function qbit_registers(code::AbstractString)
    # Split code string by lines and remove comments and preprocessor instructions.
    lines = split(code, "\n")
    filter!(s -> !(occursin("//", s) || occursin("OPENQASM 2.0", s)), lines)

    # Keep only lines containing the keyword "qreg", used to declare qbit variables.
    qreg_lines = filter(s -> occursin("qreg", s), lines)

    # Parse each line to get the register name and length: each line is of the form
    # `qreg name[n]` where `name` is an alphabetical string and `n` a positive integer.
    registers = match.(r"(?<name>[a-zA-Z0-9_]*)\[(?<length>\d+)\]", qreg_lines)
    #                              ^                   ^
    #                        register name      register length
    names = [String(r["name"]) for r in registers]
    lengths = [parse(Int, r["length"]) for r in registers]
    return names, lengths
end

"""
    qbit_sites(code::AbstractString, st::AbstractString)

Return the ITensor site indices, of SiteType `st`, associated to the quantum registers
defined in the given code.

# Example

```julia-repl
julia> code = "qreg a[3];\nqreg b[2];";

julia> qbit_sites(code, "Qubit")
5-element Vector{ITensors.Index{Int64}}:
 (dim=2|id=187|"Qubit,Site,a,n=1")
 (dim=2|id=539|"Qubit,Site,a,n=2")
 (dim=2|id=321|"Qubit,Site,a,n=3")
 (dim=2|id=981|"Qubit,Site,b,n=1")
 (dim=2|id=596|"Qubit,Site,b,n=2")
```
"""
function qbit_sites(code::AbstractString, st::AbstractString)
    # Example: with the instructions
    #   qreg a[3];
    #   qreg b[2];
    # Qiskit creates the sites a[0] a[1] a[2] b[0] b[1].
    # When `qbit_registers` is called on this code, it returns
    #   (["a", "b"], [3, 2])
    # so zip(qbit_registers(code)...) iterates over
    #   ("a", 3) ("b", 2)
    # so we end up with
    #   [siteinds(st, 3; addtags="a"); siteinds(st, 2; addtags="b")]
    return [
        [siteinds(st, n; addtags=id) for (id, n) in zip(qbit_registers(code)...)]...;
    ]
end

"""
    qbit_map(code::AbstractString)

Return a map that associates the quantum bit as written in the given
code to the relative ITensor site.

# Example

```julia-repl
julia> code = "qreg a[3];\nqreg b[2];";

julia> qbit_map(code)
Dict{String, Int64} with 5 entries:
  "a[0]" => 1
  "a[1]" => 2
  "a[2]" => 3
  "b[0]" => 4
  "b[1]" => 5
```
"""
function qbit_map(code::AbstractString)
    # Example: with the instructions
    #   qreg a[3];
    #   qreg b[2];
    # Qiskit creates the sites a[0] a[1] a[2] b[0] b[1].
    # When `qbit_registers` is called on this code, it returns
    #   (["a", "b"], [3, 2])
    #
    # This function creates a dictionary mapping Qiskit's qbit labels to a series of
    # integers, which will be the positions of each qbit in the ITensors' vector of
    # Index objects:
    #   a[0] --> 1
    #   a[1] --> 2
    #   a[2] --> 3
    #   b[0] --> 4
    #   b[1] --> 5
    names, lengths = qbit_registers(code)
    ids = []
    for (j, name) in enumerate(names)
        for k in 0:(lengths[j] - 1)
            push!(ids, "$name[$k]")
        end
    end
    # Now the map is given by the association ids[k] --> k
    return Dict(ids[k] => k for k in 1:length(ids))
end

"""
    interpret(sites::Vector{<:Index}, sitemap, line::AbstractString)

Parse `line` and return an ITensor operator derived from the OpenQASM expression, based on
the site indices in `sites`.
The `sitemap` argument is a dictionary-like map that assigns to each qbit in the OpenQASM
circuit the corresponding site number in the ITensor system. See [`qbit_map`](@ref) for
an example of such a map.

# Example

```julia-repl
julia> code = "qreg a[2];";

julia> qm = TEM.qbit_map(code);

julia> qs = TEM.qbit_sites(code, "Qubit");

julia> op = interpret(qs, qm, "cx a[0], a[1]")
ITensor ord=4 (dim=2|id=953|"Qubit,Site,a,n=2")' (dim=2|id=463|"Qubit,Site,a,n=1")' (dim=2|id=953|"Qubit,Site,a,n=2") (dim=2|id=463|"Qubit,Site,a,n=1")
NDTensors.Dense{Float64, Vector{Float64}}
```
"""
function interpret(sites::Vector{<:Index}, sitemap, line::AbstractString)
    tokens = match(r"(?<command>[^\s]*)\s(?<qbits>.*)", line)
    # Split the string at the first white space; the first part is assigned to "command",
    # the rest to "qbits". It matches
    # - [^\s]* : any sequence of non-whitespace characters
    # - .* : anything

    # Now tokens["command"] should be either a keyword or a gate name.
    if tokens["command"] == "if"
        @error "if block not (yet) implemented"
    elseif tokens["command"] == "reset"
        @error "reset instruction not (yet) implemented"
    elseif tokens["command"] == "measure"
        @error "measure instruction not (yet) implemented"
    elseif tokens["command"] == "barrier"
        @error "barrier instruction not (yet) implemented"
    elseif tokens["command"] == "qreg"
        # do nothing; we should already have dealt with registers before
    elseif tokens["command"] == "creg"
        # do nothing
    else  # it is a gate name
        gate_tokens = match(r"(?<name>[^\s()]+)(\((?<body>(?>[^()]++|(?2))*)\))", line)
        #= (from https://stackoverflow.com/a/38430207/4160978)
        Group "name" will contain the function name and "body" group will hold what is
        inside the matching parentheses.
        We need to add both ( and ) to the negated character class (?<funcion>[^\s()]+)
        because in case we have sample(decimal(3*3)) this group will grab the substring up
        to the ) (sample(decimal). Thus, we need to exclude both ( and ).

        The (\((?<body>(?>[^()]++|(?2))*)\)) part is a capture group (with ID=2) that can be
        recursed (i.e. expanded many times) with a subroutine call (?2).

        It matches
        - \( : an open round bracket
        - (?<body>(?>[^()]++|(?2))*) : Group "body" that matches zero or more sequences of:
        - [^()]++ : 1+ characters other than ( and ) or
        - (?2) : the whole \((?<body>(?>[^()]++|(?2))*)\) subpattern
        - \) : a closing parenthesis
        The (?2) subroutine call necessity (as compared to recursion with (?R)) is dictated
        by the fact that we need to repeat/recurse a part of the pattern.
        =#

        if isnothing(gate_tokens)
            # Then there were no parentheses in the string, so the gate has no parameters.
            # We redo the regex match with a simpler pattern.
            gate_tokens = match(r"(?<name>[^\s()]+)", line)
            args = []
        else
            args = eval(Meta.parse(gate_tokens[2]))
            # This will be a string contaning the parameters of the gate (angles, etc.)
        end
        gatename = gate_tokens["name"]
        # This is name of the function associated to the gate, as in qelib1_gates.jl

        qbit_list = strip.(split(tokens["qbits"], ","))
        qbit_inds = [sitemap[n] for n in qbit_list]
        qbit_siteinds = [sites[n] for n in qbit_inds]
        # `qbit_inds` contains the indices, of the `sites` list, that we need to consider.
        # `qbit_siteinds` contains the actual Index objects within `sites`.

        # Check that the arity of the gate is equal to the number or given qbits
        n_qbits = length(qbit_inds)
        gate_arity = arity(gatename)
        if (n_qbits != gate_arity)
            @error "Number of qbits ($n_qbits) ≠ arity of the gate ($gate_arity)"
        end

        if !isempty(args)
            @info "Trying to call gate $gatename with site indices" *
                "$(join(qbit_siteinds, ", ")) and additional " *
                "arguments $(join(args, ", "))"
        else
            @info "Trying to call gate $gatename with site indices " *
                "$(join(qbit_siteinds, ", ")) and no arguments"
        end

        # Get site type from `sites`
        commontags_s = ITensors.commontags(sites...)
        common_stypes = ITensors._sitetypes(commontags_s)

        for st in common_stypes
            g = gate(GateName(gatename), st, qbit_siteinds..., args...)
            !isnothing(g) && return g
        end
    end
    return nothing
end

#=
From OpenQASMexer
-----------------
'gate': [
    ('[unitary\\d+]', Token.Keyword.Type, '#push'),
    ('p\\d+', Token.Text, '#push')
],
'if_keywords': [
    ('[a-zA-Z0-9_]*', Token.Literal.String, '#pop'),
    ('\\d+', Token.Literal.Number, '#push'),
    ('.*\\(', Token.Text, 'params')
],
'index': [
    ('\\d+', Token.Literal.Number, '#pop')
],
'keywords': [
    ('\\s*("([^"]|"")*")', Token.Literal.String, '#push'),
    ('\\d+', Token.Literal.Number, '#push'),
    ('.*\\(', Token.Text, 'params')
],
'params': [
    ('[a-zA-Z_][a-zA-Z0-9_]*', Token.Text, '#push'),
    ('\\d+', Token.Literal.Number, '#push'),
    ('(\\d+\\.\\d*|\\d*\\.\\d+)([eEf][+-]?[0-9]+)?', Token.Literal.Number, '#push'),
    ('\\)', Token.Text)
],
'root': [
    ('\\n', Token.Text),
    ('[^\\S\\n]+', Token.Text),
    ('//\\n', Token.Comment),
    ('//.*?$', Token.Comment.Single),
    ('(OPENQASM|include)\\b', Token.Keyword.Reserved, 'keywords'),
    ('(qreg|creg)\\b', Token.Keyword.Declaration),
    ('(if)\\b', Token.Keyword.Reserved, 'if_keywords'),
    ('(pi)\\b', Token.Name.Constant),
    ('(barrier|measure|reset)\\b', Token.Name.Builtin, 'params'),
    ('(id|cx|x|y|z|s|sdg|h|t|tdg|ccx|c3x|c4x|c3sqrtx|rx|ry|rz|cz|cy|ch|swap|cswap|crx|cry|crz|cu1|cu3|rxx|rzz|rccx|rc3x|u1|u2|u3)\\b', Token.Keyword.Type, 'params'),
    ('[unitary\\d+]', Token.Keyword.Type),
    ('(gate)\\b', Token.Name.Function, 'gate'),
    ('[a-zA-Z_][a-zA-Z0-9_]*', Token.Text, 'index')
]
 =#

"""
    gates(code::AbstractString, st::AbstractString)

Return a pair `(s, ops)` where `s` is a list of ITensor indices of SiteType `st`, one for
each qbit declared in the given code, and `ops` is a list of ITensor operators,
representing each gate, one by one, as in the given code.

# Example

```julia-repl
julia> code = "OPENQASM 2.0;\nqreg a[2];\ncx a[0], a[1];";

julia> sites, gates = TEM.gates(code, "Qubit");

julia> sites
2-element Vector{ITensors.Index{Int64}}:
 (dim=2|id=680|"Qubit,Site,a,n=1")
 (dim=2|id=758|"Qubit,Site,a,n=2")

julia> gates[1]
ITensor ord=4 (dim=2|id=758|"Qubit,Site,a,n=2")' (dim=2|id=680|"Qubit,Site,a,n=1")' (dim=2|id=758|"Qubit,Site,a,n=2") (dim=2|id=680|"Qubit,Site,a,n=1")
NDTensors.Dense{Float64, Vector{Float64}}
```
"""
function gates(code::AbstractString, st::AbstractString)
    lines = splitlines(stripcomments(code))
    removepreprocessor!(lines)

    sites = qbit_sites(code, st)

    gates = [interpret(sites, qbit_map(code), line) for line in lines]
    filter!(!isnothing, gates)  # Remove lines unrecognized by `interpret`

    # Since `interpret` may return `nothing`, the `gates` list is actually a
    # Vector{Union{Nothing, ITensor}} object (`filter!` does not change the type), so
    # convert it (safely) to a Vector{ITensor} first.
    return sites, convert(Vector{ITensor}, gates)
end

"""
    arity(gatename::AbstractString)

Return the arity, i.e. the number of accepted/required arguments, of the gate `gatename`.
"""
function arity(gatename::AbstractString)
    arities = Dict(
        "id" => 1,
        "u1" => 1,
        "u2" => 1,
        "u3" => 1,
        "u" => 1,
        "cx" => 2,
        "x" => 1,
        "y" => 1,
        "z" => 1,
        "s" => 1,
        "p" => 1,
        "cp" => 2,
        "sdg" => 1,
        "h" => 1,
        "t" => 1,
        "tdg" => 1,
        "ccx" => 3,
        "c3x" => 4,
        "c4x" => 5,
        "rx" => 1,
        "ry" => 1,
        "rz" => 1,
        "cy" => 2,
        "cz" => 2,
        "ch" => 2,
        "swap" => 2,
        "cswap" => 3,
        "crx" => 2,
        "cry" => 2,
        "crz" => 2,
        "cu1" => 2,
        "cu3" => 2,
        "c3sqrtx" => missing,
        "rxx" => missing,
        "rzz" => missing,
        "rccx" => missing,
        "rc3x" => missing,
    )
    return arities[gatename]
end
