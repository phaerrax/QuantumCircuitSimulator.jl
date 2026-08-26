using ITensors.SiteTypes: _sitetypes, commontags
using Suppressor: @capture_err

"""
    qbitsites(code::OpenQASM.Types.MainProgram, st::AbstractString)

Return the ITensor site indices, of SiteType `st`, associated to the quantum registers
defined in the given code.

# Example

```julia-repl
julia> code = OpenQASM.parse("OPENQASM 2.0;\nqreg a[3];\nqreg b[2];");

julia> qbitsites(code, "Qubit")
5-element Vector{ITensors.Index{Int64}}:
 (dim=2|id=187|"Qubit,Site,a[0]")
 (dim=2|id=539|"Qubit,Site,a[1]")
 (dim=2|id=321|"Qubit,Site,a[2]")
 (dim=2|id=981|"Qubit,Site,b[0]")
 (dim=2|id=596|"Qubit,Site,b[2]")
```
"""
function qbitsites(code::OpenQASM.Types.MainProgram, st::AbstractString)
    registers = filter(d -> d isa OpenQASM.Types.RegDecl, code.prog)
    filter!(r -> r.type.str == "qreg", registers)  # Exclude classical registers

    sites = ITensors.Index[]
    for reg in registers
        regname = reg.name.str
        reglength = parse(Int, reg.size.str)
        append!(sites, [siteind(st; addtags=regname * "[$n]") for n in 0:(reglength - 1)])
    end
    return sites
end

function qasmstring(args::Tuple)
    stringparts = qasmstring.(args)  # Recurse into arg
    return *(stringparts...)
end

function qasmstring(args::Vector{Any})
    stringparts = qasmstring.(args)  # Recurse into arg
    return *(stringparts...)
end

function qasmstring(arg::OpenQASM.Types.Call)
    return string(arg)
end

function qasmstring(arg::RBNF.Token)
    return arg.str
end

function qasmstring(arg::OpenQASM.Types.Neg)
    return "-" * arg.val.str
end

function parsegate(sites::Vector{<:Index}, instr::OpenQASM.Types.Instruction)
    # OpenQASM.Types.Instruction fields:
    # * name: gate name
    # * cargs: numerical parameters
    # * qargs: qbits the gate acts on

    qbit_list = string.(instr.qargs)  # This will be something like ["q[1]", "q[4]"]
    qbit_inds = [getfirst(i -> hastags(i, name), sites) for name in qbit_list]
    # `qbit_inds` contains the Index objects in the `sites` list that we need to consider.
    # delete this...
    # It does _not_ contain Index objects! Altough we could write this function in an
    # equivalent way by using `getfirst` instead of `findfirst`, which would give us the
    # actual elements of the `sites` list corresponding to the indices from `findfirst`
    # (i.e. the Index objects) and then 

    parameters = @. eval(Meta.parse(qasmstring(instr.cargs)))
    if isempty(parameters)
        g = invokelatest(gate, instr.name, qbit_inds...)
    else
        g = invokelatest(gate, instr.name, qbit_inds...; cargs=parameters)
    end
    # Why `invokelatest`? The `parsegate` function is called repeatedly within `gates`; if
    # the currently parsed line is a gate declaration, then a new method is added to the
    # `gate` function (with a new GateName). However this triggers the so-called "world
    # age problem", printing the warning "method too new to be called from this world
    # context". Basically, the just-declared new method is "too new" for Julia's JIT
    # compiler to be taken into consideration. As a result, the new method basically does
    # not exist, so any following gate instruction calling that gate fails.
    # The `invokelatest` function (from Julia's Base library) allows us to circumvent this
    # issue by telling the compiler to look up the latest definition.
    return g
end

compose_txt(a, b) = "compose($a, $b)"
compose_txt(a, b...) = compose_txt(a, compose_txt(b...))

function compose(a::ITensor, b::ITensor; apply_dag::Bool=false)
    if isempty(commoninds(a, b; plev=0))
        return a * b
    else
        return apply(a, b; apply_dag=apply_dag)
    end
end

function sitetypename(::SiteType{T}) where {T}
    return string(T)
end

"""
    definition(declaration::OpenQASM.Types.Gate, st::SiteType)

Return a string containing Julia code necessary to add a `gate` method with the given name
and the site type `st` using the instructions contained in the OpenQASM `declaration`.

# Example

```julia-repl
julia> str = "OPENQASM 2.0;
qreg q[2];
gate test(a, b, c) q0, q1 {
  h q1;
  cx q0, q1;
  rz(a) q1;
  u2(b, c) q0;
  y q1;
}"

julia> g = OpenQASM.parse(str);

julia> print(QuantumCircuitSimulator.definition(g.prog[2], SiteType("Qubit")))
function QuantumCircuitSimulator.gate(::GateName"test", ::SiteType"Qubit", q0::Index, q1::Index; cargs)
a::Real = cargs[1]
b::Real = cargs[2]
c::Real = cargs[3]
compose(gate("y", q1), compose(gate("u2", q0; cargs=(b,c)), compose(gate("rz", q1; cargs=(a)), compose(gate("cx", q0, q1), gate("h", q1)))))
end
```
"""
function definition(gate::OpenQASM.Types.Gate, st::SiteType)
    # Unfortunately the `gate` method definitions require a SiteType argument, so even
    # with the "cargs as kwargs" syntax we cannot avoid passing it.
    gatename = gate.decl.name.str  # Gate name
    qbits = [qbit.str for qbit in gate.decl.qargs]  # Gate qubit arguments
    params = [p.str for p in gate.decl.cargs]  # Gate parameter names

    index_list = ["$q::Index" for q in qbits]

    fn_signature_decl = "function QuantumCircuitSimulator.gate("
    fn_signature_parts = [
        "::GateName\"" * gatename * "\""
        "::SiteType\"" * sitetypename(st) * "\""
        index_list
    ]
    fn_signature = if isempty(params)
        fn_signature_decl * join(fn_signature_parts, ", ") * ")"
    else
        fn_signature_decl * join(fn_signature_parts, ", ") * "; cargs)"
    end

    fn_cargs_declaration = ["$p::Real = cargs[$i]" for (i, p) in enumerate(params)]

    fn_body = []
    for instr in gate.body
        name = "\"" * instr.name * "\""
        indices = [qa.name.str for qa in instr.qargs]
        parameters = [qasmstring(ca) for ca in instr.cargs]
        if isempty(parameters)
            push!(fn_body, "gate(" * join([name; indices], ", ") * ")")
        else
            push!(
                fn_body,
                "gate(" *
                join([name; indices], ", ") *
                "; cargs=(" *
                join(parameters, ",") *
                "))",
            )
        end
    end

    # We transform the body so that the gates are correctly multiplied. We use the
    # `apply` function, i.e. `[a, b, c]` is turned into `apply(a, apply(b, c))`.
    if length(fn_body) == 1  # TODO: what if length(fn_body) == 0?
        str = join([fn_signature; fn_cargs_declaration; fn_body; "end"], "\n")
    else
        str = join(
            [fn_signature; fn_cargs_declaration; compose_txt(reverse(fn_body)...); "end"],
            "\n",
        )
    end
    return str
end

function eval_gate_definition(code::AbstractString; warn_on_gate_redefinition)
    expr = Meta.parse(code)
    # When `expr` is finally evaluated, it can trigger warnings such as
    #
    # WARNING: Method definition gate(QuantumCircuitSimulator.GateName{:rzx},
    # ITensors.SiteTypes.SiteType{ITensors.SmallStrings.SmallString(data=StaticArraysCore.SArray{Tuple{16},
    # UInt16, 1, 16}(data=(0x0051, 0x0075, 0x0062, 0x0069, 0x0074, 0x0000, 0x0000, 0x0000,
    # 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000)))}, ITensors.Index{T}
    # where T, ITensors.Index{T} where T) in module QuantumCircuitSimulator at none:1
    # overwritten on the same line (check for duplicate calls to `include`).
    #
    # due to the fact that the current gate definition ends up overwriting an existing one.
    # The warning is printed even if the old and the new methods are identical, i.e. it is
    # shown whenever the method signatures are identical.
    # This is usually not a problem, because the input OpenQASM file should be taken as the
    # ultimate source of truth for gate definitions, so even when a gate from the input file
    # ends up overwriting one in the package library, there is a good reason for it: the
    # circuit simulation should follow the definition given in the OpenQASM file.
    #
    # For this reason, these warnings are mostly annoying and can be suppressed. We keep
    # however the option to show them anyway, it might be useful for some purposes.

    if warn_on_gate_redefinition
        # * Option A: show the warnings anyway. This is standard Julia execution.
        return eval(expr)
    else
        # * Option B: suppress the warnings.
        # Since they are emitted by Julia's compiler, and printed straight to stderr, we
        # capture the stream temporarily. However, this can also inadvertently suppress
        # anything else that is passed to stderr while it is captured. For this reason,
        # in order not to lose that output, we re-print it later, filtering out the method
        # overwrite warnings.
        # The overwrite warning has a fixed shape: it always starts with "WARNING: Method
        # definition" and always ends with "(check for duplicate calls to `include`).", so
        # we can parse the whole message as one block with a regex, then strip away the
        # match, leaving everything else untouched.
        msgs = @capture_err eval(expr)
        filtered_msgs = replace(
            msgs,
            r"WARNING: Method definition.*?\(check for duplicate calls to `include`\)\.\n?"s => "",
            # - the `s` flag at the end lets the dot (.) match across newlines as well ---
            #   useful in case the warning itself wraps.
            # - the `.*?` sequence between the beginning and end snippets ensures that each
            #   match stops at the nearest closing phrase (it's a non greedy quantifier).
        )
        print(stderr, filtered_msgs)
    end

    return nothing
end

function _parsegates(
    code::OpenQASM.Types.MainProgram,
    sites::Vector{<:Index},
    st::AbstractString;
    warn_on_gate_redefinition,
)
    gatelist = ITensor[]
    for line in code.prog
        if line isa OpenQASM.Types.Instruction
            push!(gatelist, parsegate(sites, line))
        elseif line isa OpenQASM.Types.Gate
            new_gate_method = definition(line, SiteType(st))
            eval_gate_definition(new_gate_method; warn_on_gate_redefinition)
        end
    end
    return gatelist
end

"""
    gates(code::OpenQASM.Types.MainProgram, st::AbstractString; warn_on_gate_redefinition)

Create a list of gates (ITensor operators) as parsed from the given OpenQASM `code`,
returning a tuple `(s, ops)` where:

* `s` is a vector of ITensor indices of SiteType `st` (an index for each qbit declared in
  `code`),
* `ops` is a list of ITensor operators, one for each gate, in the appearance order in `code`.

# Example

```julia-repl
julia> code = OpenQASM.parse("OPENQASM 2.0;\nqreg a[2];\ncx a[0], a[1];");

julia> s, g = gates(code, "Qubit"; warn_on_gate_redefinition=false);

julia> s
2-element Vector{ITensors.Index{Int64}}:
 (dim=2|id=680|"Qubit,Site,a,n=1")
 (dim=2|id=758|"Qubit,Site,a,n=2")

julia> g[1]
ITensor ord=4 (dim=2|id=758|"Qubit,Site,a,n=2")' (dim=2|id=680|"Qubit,Site,a,n=1")' (dim=2|id=758|"Qubit,Site,a,n=2") (dim=2|id=680|"Qubit,Site,a,n=1")
NDTensors.Dense{Float64, Vector{Float64}} ```
"""
function gates(
    code::OpenQASM.Types.MainProgram, st::AbstractString; warn_on_gate_redefinition
)
    sites = qbitsites(code, st)
    gatelist = _parsegates(code, sites, st; warn_on_gate_redefinition)
    return sites, gatelist
end

"""
    gates(code::OpenQASM.Types.MainProgram, sites::Vector{<:Index})

Return a list of gates (ITensor operators) as parsed from the given OpenQASM `code`,
building the operators on the already existing Index objects in `sites`.
The given `sites` are checked to ensure they are compatible to the ones that
would be generated from `code`.
"""
function gates(
    code::OpenQASM.Types.MainProgram, sites::Vector{<:Index}; warn_on_gate_redefinition
)
    # Determine whether _all_ indices in `sites` are of the Qubit or vQubit site type.
    # The site type of these indices will determine the site type of the gates as well.
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    if "Qubit" in sitetypename.(common_stypes)
        st = "Qubit"
    elseif "vQubit" in sitetypename.(common_stypes)
        st = "vQubit"
    else
        error("unsupported site type: only \"Qubit\" and \"vQubit\" are allowed.")
    end

    # Read how many qbits are declared by the input OpenQASM program, and check that there
    # are as many site indices in `sites`.
    newsites = qbitsites(code, st)
    if length(newsites) != length(sites)
        error("input OpenQASM program declares $(length(newsites)) qbits, but " *
              "only $(length(sites)) were provided to this function.")
    end

    return _parsegates(code, sites, st; warn_on_gate_redefinition)
end
