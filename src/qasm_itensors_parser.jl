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
    qbit_siteinds = [sites[n] for n in qbit_inds]
    # `qbit_inds` contains the indices, of the `sites` list, that we need to consider.
    # `qbit_siteinds` contains the actual Index objects within `sites`.

    # Get site type from `sites`
    commontags_s = ITensors.commontags(sites...)
    common_stypes = ITensors._sitetypes(commontags_s)

    parameters = @. eval(Meta.parse(qasmstring(instr.cargs)))
    for st in common_stypes
        g = gate(GateName(instr.name), st, qbit_siteinds..., parameters...)
        !isnothing(g) && return g
    end
    return nothing
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

    gates =
        parsegate.(
            Ref(sites),
            Ref(qmap),
            filter(line -> line isa OpenQASM.Types.Instruction, code.prog),
        )

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

    layers = MPO[]
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
            push!(layers, MPO(currentlayer))
            currentlayer = [pop!(gatestack)]
        else
            # If not: add the gate to the current layer and go on.
            push!(currentlayer, pop!(gatestack))
        end
    end

    return layers
end