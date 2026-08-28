const GateLayer = Vector{ITensor}

function ITensorMPS.MPO(sites::Vector{<:ITensors.Index}, gl::GateLayer; kwargs...)
    # TODO this is a lazy implementation, we could do something better and do the
    # decomposition here directly, instead of multiplying.
    m = MPO(sites, "Id")
    for g in gl
        m = apply(g, m; kwargs...)
    end

    return m
end

"""
    QuantumCircuit

A `QuantumCircuit` is an aggregate type that stores the qbit array of the circuit and the
sequence of gates, or instructions, arranged in layers.
"""
struct QuantumCircuit
    sites::Vector{<:ITensors.Index}
    instructions::Vector{GateLayer}
    function QuantumCircuit(
        sites::Vector{<:ITensors.Index}, instructions::Vector{GateLayer}
    )
        # Sanity check: each gate layer in `instructions` contains ITensors defined on
        # `sites`, and no other indices.
        for (n, gl) in enumerate(instructions)
            bad_qbits = setdiff(occupiedsites(gl), sites)
            if !isempty(bad_qbits)
                error(
                    "gates and sites do not match: layer $n contains gate(s) acting on " *
                    "$(join(_qbittag.(bad_qbits), ", ")), which are not among the " *
                    "circuit sites.",
                )
            end
        end
        return new(sites, instructions)
    end
end

ITensorMPS.siteinds(circ::QuantumCircuit) = circ.sites
instructions(circ::QuantumCircuit) = circ.instructions

"""
    QuantumCircuit(sites)

Return an empty circuit, in which the qbit array is defined but which contains no gates. 
"""
function QuantumCircuit(sites)
    return QuantumCircuit(sites, GateLayer[])
end

depth(circ::QuantumCircuit) = length(circ.instructions)
nqbits(circ::QuantumCircuit) = length(circ.sites)

# Iteration utilities
# -------------------
# Iterating over the circuit == iterating over its instructions (i.e. its layers)
Base.length(circ::QuantumCircuit) = depth(circ)  # required for enumerating
Base.eachindex(circ::QuantumCircuit) = eachindex(circ.instructions)
Base.getindex(circ::QuantumCircuit, idx) = getindex(circ.instructions, idx)
Base.lastindex(circ::QuantumCircuit) = lastindex(circ.instructions)
Base.iterate(circ::QuantumCircuit) = iterate(circ.instructions)
Base.iterate(circ::QuantumCircuit, state) = iterate(circ.instructions, state)
Base.push!(circ::QuantumCircuit, instr) = push!(circ.instructions, instr)

"""
    occupiedsites(gl::GateLayer)

Return the indices within `sites` on which the gates in the `gl` layer act.
"""
function occupiedsites(gl::GateLayer)
    return collect(Iterators.flatten(inds.(gl; plev=0)))
end

"""
    instructionsites(instr::OpenQASM.Types.Instruction, sites::Vector{<:Index})
    instructionsites(instr::OpenQASM.Types.Barrier, sites::Vector{<:Index})

Return the indices within `sites` on which the instruction or barrier `instr` acts.
"""
function instructionsites(
    instr::Union{OpenQASM.Types.Instruction,OpenQASM.Types.Barrier},
    sites::Vector{<:ITensors.Index},
)
    sitelist = ITensors.Index[]
    for qbit in string.(instr.qargs)
        append!(sitelist, filter(hastags(qbit), sites))
    end
    return unique(sitelist)  # Probably `unique` is unnecessary here, but just in case...
end

"""
    freesites(sites::Vector{<:ITensor.Index}, gl::GateLayer)

Return a list of the sites in `sites` that are not acted upon by any gate from the `gl`
layer.
"""
function freesites(sites::Vector{<:ITensors.Index}, gl::GateLayer)
    return setdiff(sites, occupiedsites(gl))
end

"""
    freesites(circ::QuantumCircuit, gl::GateLayer)

Return a list of the sites in the circuit `circ` with no gates from the `gl` layer acting
upon them.
"""
freesites(circ::QuantumCircuit, gl::GateLayer) = freesites(siteinds(circ), gl)

"""
    lastblockinglayer(circ::QuantumCircuit, sites)

Return the index of the layer that something (a gate, or a barrier) acting on `sites` cannot
slide past: it may be placed in any layer after this one, but not in this one or any earlier
one. Return `0` if it could be placed in the very first layer (in particular, if `circ` is
empty).
"""
function lastblockinglayer(circ::QuantumCircuit, sites)
    # We need to decide which layer this belongs into. It will need to be placed in the
    # "earliest" layer such that it itself _and all later layers_ have free spots in each of
    # `sites`. (This should mimic Qiskit's "asap" layering algorithm.) In other words, it
    # can share a layer with earlier gates only if none of them touch `sites`, and clearly
    # it can't be placed before an earlier gate that touches one of its qbits.  So, starting
    # from the newest layer and walking backward, it can keep sliding earlier as long as
    # `sites` are completely free in every layer it passes through, and it stops at the
    # first (i.e., most recent) layer where one of its qbits is already occupied or at a
    # barrier.

    idx = 0
    @debug "Sites: " * join(_qbittag.(sites), " ")
    # We start from assuming `idx == 0` and then in the following loop we scan the already
    # existing layers to find out possible obstructions. If no obstruction is found, then
    # the loop never breaks (all existing layers are free for `sites`, or `circ` is empty)
    # and `idx` stays zero, signalling that it can be placed in the very first layer.

    for (j, gl) in Iterators.reverse(enumerate(circ))
        # To check whether it fits in each of the latest group of layers. Start from the last
        # of the already formed layers and go backwards, finding the latest one where it
        # doesn't fit: the _next_ one will be the layer it needs to go in.
        actual_layer_idx = j
        @debug "Checking layer $actual_layer_idx. Free spots: " *
            join(_qbittag.(freesites(circ, gl)), " ")
        if !issubset(sites, freesites(circ, gl))
            # Then layer `j` already uses one of the qbits, so the layer is the most recent
            # one that blocks it; consequently, it must go in layer `j+1` (the layer right
            # after the blocker).
            @debug "No space available on this layer."
            idx = j
            break
        end
        # else: the qubits are all free in this layer, so we keep looking further back.
    end

    return idx
end

"""
    split_by(f, a::Vector; keeptrue=true)

Split the vector `a` in subvectors on wherever `f` is `true`. If `keeptrue` is true, the
elements of `a` such that `f(a)` is `true` are kept, and will be placed at the end of each
subvector.

# Source

https://stackoverflow.com/a/63739895/4160978
"""
function split_by(f, a::Vector; keeptrue=true)
    result = Vector{typeof(view(a, 1:0))}()
    l = firstindex(a)
    r = firstindex(a)
    while r <= lastindex(a)
        if f(a[r])
            if keeptrue
                push!(result, @view(a[l:r]))
            else
                push!(result, @view(a[l:(r - 1)]))
            end
            l = r + 1
        end
        r += 1
    end
    push!(result, @view(a[l:end]))

    return result
end

function _qbittag(s)
    taglist = string.(collect(tags(s)))
    idx = findfirst(t -> contains(t, "[") && contains(t, "]"), taglist)
    return taglist[idx]
end

"""
    QuantumCircuit(code::OpenQASM.Types.MainProgram; kwargs...)
    QuantumCircuit(code::AbstractString; kwargs...)

Parse the circuit defined in `code` (either a string or an already parsed OpenQASM program)
as a `QuantumCircuit` object, structuring the sequence of gates in layer. The gates are
added to the layers in a way that _should_ resemble the "ASAP" policy by Qiskit.

# Keyword arguments

* `operator_picture`: if `false` (default) then the circuit is interpreted as a sequence of
  gates to be applied to a pure state; if it is `true` instead it will be a sequence of
  gates to be applied to mixed states (in the Schrödinger picture) or to observables (in the
  Heisenberg picture).
* `warn_on_gate_redefinition` if `false` (default), ignore the "Method overwritten" warnings
  that the Julia compiler might print whenever a gate definition overwrites an existing one.
"""
function QuantumCircuit(
    code::OpenQASM.Types.MainProgram;
    operator_picture=false,
    warn_on_gate_redefinition=false,
)
    st = (operator_picture ? "vQubit" : "Qubit")
    circsites = qbitsites(code, st)

    circ = QuantumCircuit(circsites)

    # Maps each site to the earliest layer index that a gate acting on it is allowed to be
    # placed after, due to a previously encountered barrier on that site. A barrier doesn't
    # add anything to `circ` itself: it just raises this floor for the qbits it acts on, so
    # that later gates on those qbits can no longer go back past it.
    barrier_floor = Dict{ITensors.Index,Int}()

    for line in code.prog
        if line isa OpenQASM.Types.Instruction
            g = parsegate(circsites, line)  # build the ITensor for this gate
            g_sites = inds(g; plev=0)  # qbits/sites the ITensors acts on

            # `g` must go in the layer right after the last one that blocks it (see
            # `lastblockinglayer`), i.e. layer `idx + 1`.
            idx = lastblockinglayer(circ, g_sites)

            # A previous barrier on one of `g`'s qbits may forbid it from landing as early
            # as `idx`: we raise `idx` to the highest barrier floor among `g_sites`, if
            # there is any.
            barrier_idx = maximum(get(barrier_floor, s, 0) for s in g_sites; init=0)
            idx = max(idx, barrier_idx)

            if idx + 1 > depth(circ)
                # This means that the first "blocking" gate is already the last one (or
                # `circ` is empty), so we need to add a new layer to the circuit.
                @debug "New layer required for gate \"$line\"."
                push!(circ, [g])
                @debug "Added new layer at the end: new circuit length is $(depth(circ))."
            elseif idx == 0
                # The loop ran until the end: `g` fits in all layers in `circ`, which means
                # that it will end up in the very first layer.
                # What if this is the first gate and the circuit is still empty?
                # If `circ` is empty, then `idx + 1` will be greater than `depth(circ)`,
                # which is zero, so the previous case will be selected. We don't need to
                # worry here whether there actually is a "first layer" which we can put
                # `g` in.
                @debug "Reached the beginning of the circuit. " *
                    "Added gate \"$line\" to layer $(idx+1)."
                push!(circ[idx + 1], g)
            else
                @debug "Added gate \"$line\" to layer $(idx+1)."
                push!(circ[idx + 1], g)
            end
        elseif line isa OpenQASM.Types.Gate
            # This is a definition of a new gate, so we call the code that creates and
            # evaluates a new `gate` method.
            new_gate_method = definition(line, SiteType(st))
            eval_gate_definition(new_gate_method; warn_on_gate_redefinition)
        elseif line isa OpenQASM.Types.Barrier
            # We encountered a barrier: get the qbits affected by the barrier, and raise the
            # respective floor. A barrier is not a gate and never occupies a layer by
            # itself, but it should still slide backwards as far as possible, like gates do.
            # We also take into account any floor a previous, overlapping barrier might have
            # already set on these qbits.
            barrier_sites = instructionsites(line, circsites)
            idx = lastblockinglayer(circ, barrier_sites)
            idx = max(idx, maximum(get(barrier_floor, s, 0) for s in barrier_sites; init=0))
            for s in barrier_sites
                barrier_floor[s] = idx
            end
        end
    end

    return circ
end

function QuantumCircuit(code::AbstractString; kwargs...)
    return QuantumCircuit(OpenQASM.parse(code); kwargs...)
end

"""
    layers_mpo(circ; progress=false)

Return a list of MPOs, one for each layer of the circuit `circ`.
Optionally display a progress bar by setting `progress` to `true`.
"""
function layers_mpo(circ::QuantumCircuit; progress=false)
    # The layers have already been formed while parsing: it's just a matter of converting
    # them into MPOs.
    mpos = Vector{MPO}(undef, depth(circ))
    if progress
        pbar = Progress(depth(circ); desc="Creating MPO from layers...")
    end
    for (i, ℓ) in enumerate(circ)
        mpos[i] = MPO(siteinds(circ), ℓ)
        progress && next!(pbar)
    end
    return mpos
end
