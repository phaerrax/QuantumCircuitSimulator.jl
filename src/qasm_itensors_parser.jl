function qbitregisters(code::OpenQASM.Types.MainProgram)
    registers = filter(d -> d isa OpenQASM.Types.RegDecl, code.prog)
    filter!(r -> r.type.str == "qreg", registers)  # Exclude classical registers
    return [(r.name.str, parse(Int, r.size.str)) for r in registers]
end

"""
    qbitsites(code::OpenQASM.Types.MainProgram, st::AbstractString)

Return the ITensor site indices, of SiteType `st`, associated to the quantum registers
defined in the given code.

# Example

```julia-repl
julia> code = OpenQASM.parse("qreg a[3];\nqreg b[2];");

julia> qbitsites(code, "Qubit")
5-element Vector{ITensors.Index{Int64}}:
 (dim=2|id=187|"Qubit,Site,a,n=1")
 (dim=2|id=539|"Qubit,Site,a,n=2")
 (dim=2|id=321|"Qubit,Site,a,n=3")
 (dim=2|id=981|"Qubit,Site,b,n=1")
 (dim=2|id=596|"Qubit,Site,b,n=2")
```
"""
function qbitsites(code::OpenQASM.Types.MainProgram, st::AbstractString)
    return [
        [siteinds(st, n; addtags=id) for (id, n) in qbitregisters(code)]...;
    ]
end

"""
    qbitmap(code::OpenQASM.Types.MainProgram)

Return a map that associates the quantum bit as written in the given
code to the relative ITensor site.

# Example

```julia-repl
julia> code = OpenQASM.parse("qreg a[3];\nqreg b[2];");

julia> qbitmap(code)
Dict{String, Int64} with 5 entries:
  "a[0]" => 1
  "a[1]" => 2
  "a[2]" => 3
  "b[0]" => 4
  "b[1]" => 5
```
"""
function qbitmap(code::OpenQASM.Types.MainProgram)
    # Example: with the instructions
    #   qreg a[3];
    #   qreg b[2];
    # Qiskit creates the sites a[0] a[1] a[2] b[0] b[1].
    # When `qbitregisters` is called on this code, it returns
    #   [("a", 3), ("b", 2)]
    #
    # This function creates a dictionary mapping Qiskit's qbit labels to a series of
    # integers, which will be the positions of each qbit in the ITensors' vector of
    # Index objects:
    #   a[0] --> 1
    #   a[1] --> 2
    #   a[2] --> 3
    #   b[0] --> 4
    #   b[1] --> 5
    registers = qbitregisters(code)
    ids = []
    for (name, len) in registers
        for k in 0:(len - 1)
            push!(ids, "$name[$k]")
        end
    end
    # Now the map is given by the association ids[k] => k
    return Dict(ids[k] => k for k in 1:length(ids))
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

function parsegate(sites::Vector{<:Index}, sitemap, instr::OpenQASM.Types.Instruction)
    # OpenQASM.Types.Instruction fields:
    # * name: gate name
    # * cargs: numerical parameters
    # * qargs: qbits the gate acts on

    qbit_list = string.(instr.qargs)  # This will be something like ["q[1]", "q[4]"]
    qbit_inds = [sitemap[n] for n in qbit_list]
    # `qbit_inds` contains the indices, of the `sites` list, that we need to consider.

    parameters = @. eval(Meta.parse(qasmstring(instr.cargs)))
    if isempty(parameters)
        g = invokelatest(gate, instr.name, sites, qbit_inds...)
    else
        g = invokelatest(gate, instr.name, sites, qbit_inds...; cargs=parameters)
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
        return apply(a, b; apply_dag)
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

julia> print(TEM.definition(g.prog[2], SiteType("Qubit")))
function TEM.gate(::GateName"test", ::SiteType"Qubit", q0::Index, q1::Index; cargs)
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

    fn_signature_decl = "function TEM.gate("
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

Return a pair `(s, ops)` where `s` is a list of ITensor indices of SiteType `st`, one for
each qbit declared in the given code, and `ops` is a list of ITensor operators,
representing each gate, one by one, as in the given code.

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
    qmap = qbitmap(code)

    gates = []
    for line in code.prog
        if line isa OpenQASM.Types.Instruction
            push!(gates, parsegate(sites, qmap, line))
        elseif line isa OpenQASM.Types.Gate
            new_gate_method = definition(line, SiteType(st))
            eval(Meta.parse(new_gate_method))
        end
    end
    return sites, gates
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
        if havecommonelements(
            collect(inds(first(gatestack))), Iterators.flatten(inds.(currentlayer))
        )
            # If yes: stop adding to the layer, save it and start a new layer with this gate
            push!(layers, currentlayer)
            currentlayer = [pop!(gatestack)]
        else
            # If not: add the gate to the current layer and go on.
            push!(currentlayer, pop!(gatestack))
        end
    end

    return layers
end
