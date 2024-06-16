using ITensors.SiteTypes: _sitetypes, commontags

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

"""
    gates(code::OpenQASM.Types.MainProgram, st::AbstractString)

Create a list of gates (ITensor operators) as parsed from the given OpenQASM `code`,
returning a triple `(s, ops)` where:

* `s` is a list of ITensor indices of SiteType `st`, one for each qbit declared in the
given code,
* `ops` is a list of ITensor operators, one for each gate, in order.

# Example

```julia-repl
julia> code = OpenQASM.parse("OPENQASM 2.0;\nqreg a[2];\ncx a[0], a[1];");

julia> s, g = gates(code, "Qubit");

julia> s
2-element Vector{ITensors.Index{Int64}}:
 (dim=2|id=680|"Qubit,Site,a,n=1")
 (dim=2|id=758|"Qubit,Site,a,n=2")

julia> g[1]
ITensor ord=4 (dim=2|id=758|"Qubit,Site,a,n=2")' (dim=2|id=680|"Qubit,Site,a,n=1")' (dim=2|id=758|"Qubit,Site,a,n=2") (dim=2|id=680|"Qubit,Site,a,n=1")
NDTensors.Dense{Float64, Vector{Float64}}
```
"""
function gates(code::OpenQASM.Types.MainProgram, st::AbstractString)
    sites = qbitsites(code, st)

    gatelist = ITensor[]
    for line in code.prog
        if line isa OpenQASM.Types.Instruction
            push!(gatelist, parsegate(sites, line))
        elseif line isa OpenQASM.Types.Gate
            new_gate_method = definition(line, SiteType(st))
            eval(Meta.parse(new_gate_method))
        end
    end
    return sites, gatelist
end

"""
    gates(code::OpenQASM.Types.MainProgram, sites::Vector{<:Index})

Return a list of gates (ITensor operators) as parsed from the given OpenQASM `code`,
building the operators on the already existing Index objects in `sites`.
The given `sites` are checked to ensure they are compatible to the ones that
would be generated from `code`.
"""
function gates(code::OpenQASM.Types.MainProgram, sites::Vector{<:Index})
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    if "Qubit" in sitetypename.(common_stypes)
        st = "Qubit"
    elseif "vQubit" in sitetypename.(common_stypes)
        st = "vQubit"
    else
        error("Unrecognized SiteType: only \"Qubit\" and \"vQubit\" are allowed.")
    end

    newsites = qbitsites(code, st)
    if length(newsites) != length(sites)
        error("Sites not compatible with input OpenQASM code.")
    end

    gates = ITensor[]
    for line in code.prog
        if line isa OpenQASM.Types.Instruction
            push!(gates, parsegate(sites, line))
        elseif line isa OpenQASM.Types.Gate
            new_gate_method = definition(line, SiteType(st))
            eval(Meta.parse(new_gate_method))
        end
    end
    return gates
end

"""
    gatelayers(gates::Vector{ITensor})

Return a list of MPOs, each one representing a layer of the circuit.
A layer is created multiplying together as many adjacent gates as possible, i.e. the
construction of the MPO ends as soon as the next gate in line acts on a qbit which is
already acted on by the gates already in the MPO.
This way, all gates in a layer commute with each other and can be executed in any order.
"""
function gatelayers(gates::Vector{ITensor})
    # We can detect this by looking at the Indices of the gate ITensors: as soon as we see a
    # repeated index, we stop grouping them and start a new layer.
    # TODO: Implement circuit barriers.
    # We put the gates in a stack, so that we can pop them one at a time. We reverse the
    # order, so that the first gate in the list is also the first in the stack.
    gatestack = Stack{ITensor}()
    foreach(g -> push!(gatestack, g), reverse(gates))

    layers = Vector{ITensor}[]
    # Will contain only gates, no noise MPOs. Each layer will be a vector of gates.

    # Create a new empty layer.
    currentlayer = ITensor[]
    while !isempty(gatestack)
        # Look up the next tensor in the gate sequence: it is `first(gatestack)`.
        # Does this gate have an Index which is already present in the current layer?
        # 1. inds(first(gatestack)) is a tuple of Index objects
        # 2. inds.(currentlayer) is a _vector_ of tuples of Index objects
        # we need both of them to be a big list of Index, not of tuples:
        # 1. `collect` transforms a Tuple into a Vector
        # 2. `flatten` transforms a Vector of Tuple{X,X,...} into a Vector{X}
        nextgate = pop!(gatestack)

        inds_next = collect(inds(nextgate))  # indices of the next gate in the queue
        inds_current = Iterators.flatten(inds.(currentlayer))  # all indices in the layer
        inds_common = intersect(inds_next, inds_current)
        if isempty(inds_common)
            # If there are no common indices: add the gate to the current layer and go on.
            push!(currentlayer, nextgate)
        else
            # If there are: stop adding to the layer, save it and start a new layer with
            # this gate.
            push!(layers, currentlayer)
            currentlayer = [nextgate]
        end
    end
    # When the gate stack is finally empty, the current layer is yet to be pushed.
    push!(layers, currentlayer)

    return layers
end
