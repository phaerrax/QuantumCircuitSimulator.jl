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
