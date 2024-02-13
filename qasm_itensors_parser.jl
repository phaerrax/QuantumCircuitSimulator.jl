using ITensors

"""
    splitlines(code::AbstractString)

Split code at semicolons, producing a list whose elements are single lines of code.
"""
function splitlines(str)
    lines = split(str, ";")  # Split at semicolons
    lines = strip.(lines)  # Remove trailing and leading whitespace (such as stray '\n's)
    return filter(!isempty, lines)  # Remove empty lines
end

function stripcomments(str)
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
    qbit_sites(code::AbstractString)

Return names and lengths of the quantum registers declared in the given
code.
"""
function qbit_registers(code::AbstractString)
    lines = split(code, "\n")
    filter!(s -> !(occursin("//", s) || occursin("OPENQASM 2.0", s)), lines)
    qreg_lines = filter(s -> occursin("qreg", s), lines)
    registers = match.(r"(?<name>[a-zA-Z0-9_]*)\[(?<length>\d+)\]", qreg_lines)
    #                              ^                   ^
    #                     variable name      register length
    names = [r["name"] for r in registers]
    lengths = [parse(Int, r["length"]) for r in registers]
    return String.(names), lengths
end

"""
    qbit_sites(code::AbstractString)

Return the ITensor site Index vector associated to the quantum registers
defined in the given code.
"""
function qbit_sites(code::AbstractString)
    return [
        [siteinds("Qubit", n; addtags=id) for (id, n) in zip(qbit_registers(code)...)]...;
    ]
end

"""
    qbit_map(code::AbstractString)

Return a map that associates the quantum bit as written in the given
code to the relative ITensor site.
"""
function qbit_map(code::AbstractString)
    # Example:
    #   qreg q[3]
    #   qreg ancilla[1]
    # results in
    #   ("q", 3)
    #   ("ancilla", 1)
    # when `qbit_registers` is called. The map should read
    # "q[0]"       --> 1
    # "q[1]"       --> 2
    # "q[2]"       --> 3
    # "ancilla[0]" --> 4
    names, lengths = qbit_registers(code)
    ids = []
    for (j, name) in enumerate(names)
        for k in 0:(lengths[j] - 1)
            push!(ids, "$name[$k]")
        end
    end
    # Now the map is given by ids[k] --> k
    return Dict(ids[k] => k for k in 1:length(ids))
end

"""
    interpret(sites::Vector{<:Index}, sitemap, line::AbstractString)

Parse `line` and return an ITensor operator derived from the OpenQASM expression, based on
the site indices in `sites`.
The `sitemap` argument is a dictionary-like map that assigns to each qbit in the OpenQASM
circuit the corresponding site number in the ITensor system. See [`qbit_map`](@ref) for
an example of such a map.
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
        end
        gate_fn = "gate_" * gate_tokens["name"]
        # This is name of the function associated to the gate, as in qelib1_gates.jl

        qbit_list = strip.(split(tokens["qbits"], ","))
        qbit_siteinds = [sitemap[n] for n in qbit_list]

        # Check that the arity of the gate is equal to the number or given qbits
        n_qbits = length(qbit_siteinds)
        gate_arity = arity(replace(gate_fn, "gate_" => ""))
        if (n_qbits != gate_arity)
            @error "Number of qbits ($n_qbits) ≠ arity of the gate ($gate_arity)"
        end

        if !isempty(args)
            @info "Trying to call $gate_fn with site indices $(join(qbit_siteinds, ", ")) " *
                "and additional arguments $(join(args, ", "))"
        else
            @info "Trying to call $gate_fn with site indices $(join(qbit_siteinds, ", ")) " *
                "and no arguments"
        end
        # This will be a string contaning the parameters of the gate (angles, etc.)
        return getfield(Main, Symbol(gate_fn))(sites, qbit_siteinds..., args...)
        # Finally, call the function given by `gate_fn` with all the needed arguments
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
    gates(code::AbstractString)

Return a pair `(s, ops)` where `s` is a list of ITensor "Qubit" indices, one for each qbit
declared in the given code, and `ops` is a list of ITensor operators, representing each
gate, one by one, as in the given code.
"""
function gates(code::AbstractString)
    lines = splitlines(stripcomments(code))
    removepreprocessor!(lines)

    s = qbit_sites(code)

    gates = interpret.(Ref(s), Ref(qbit_map(code)), lines)
    filter!(!isnothing, gates)  # Remove lines unrecognized by `interpret`
    # Since `interpret` may return `nothing`, the `gates` list is actually a
    # Vector{Union{Nothing, ITensor}} object (`filter!` does not change the type).

    return s, convert(Vector{ITensor}, gates)
end
