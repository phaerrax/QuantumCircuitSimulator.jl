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
end

ITensorMPS.siteinds(circ::QuantumCircuit) = circ.sites
instructions(circ::QuantumCircuit) = circ.instructions

"""
    quantumcircuit(sites)

Return an empty circuit, in which the qbit array is defined but which contains no gates. 
"""
function quantumcircuit(sites)
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
    findsites(instr::OpenQASM.Types.Instruction, sites::Vector{<:Index})

Return the list of indices within `sites` on which the instruction `instr` acts.
"""
function findsites(instr::OpenQASM.Types.Instruction, sites::Vector{<:ITensors.Index})
    sitelist = Int[]
    for qbit in string.(instr.qargs)
        append!(sitelist, findall(hastags(qbit), s))
    end
    return sitelist
end

"""
    freesites(sites::Vector{<:ITensor.Index}, gl::GateLayer)

Return a list of the sites in `sites` that are not acted upon by any gate from the `gl`
layer.
"""
function freesites(sites::Vector{<:ITensors.Index}, gl::GateLayer)
    return setdiff(sites, collect(Iterators.flatten(inds.(gl; plev=0))))
end

"""
    freesites(circ::QuantumCircuit, gl::GateLayer)

Return a list of the sites in the circuit `circ` with no gates from the `gl` layer acting
upon them.
"""
freesites(circ::QuantumCircuit, gl::GateLayer) = freesites(siteinds(circ), gl)

"""
    freesites(circ::QuantumCircuit)

Return a dictionary of the sites in each layer of the circuit `circ` with no gates, in the
circuit itself, acting upon them.
"""
function freesites(circ::QuantumCircuit)
    return Dict(i => freesites(circ, gl) for (i, gl) in enumerate(circ))
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
    quantumcircuit(code::OpenQASM.Types.MainProgram; operator_picture=false)

Parse the circuit defined in `code` as a `QuantumCircuit` object, structuring the sequence
of gates in layer. The gates are added to the layers in a way that _should_ resemble
the "ASAP" policy by Qiskit.

If `operator_picture` is `false` then the circuit is interpreted as a sequence of gates to
be applied to a pure state; if it is `true` instead it will be a sequence of gates to be
applied to mixed states (in the Schrödinger picture) or to observables (in the Heisenberg
picture).
"""
function quantumcircuit(code::OpenQASM.Types.MainProgram; operator_picture=false)
    st = (operator_picture ? "vQubit" : "Qubit")
    circsites = qbitsites(code, st)

    circ = quantumcircuit(circsites)

    for line in code.prog
        if line isa OpenQASM.Types.Instruction
            # Parse the gate instructions
            g = parsegate(circsites, line)

            g_sites = inds(g; plev=0)  # sites on which the gate is defined
            # The gate will need to be placed in the "earliest" layer such that
            # itself _and all later layers_ have free spots in each site in `g_sites`.
            # This should mimic Qiskit's "asap" layering algorithm.

            # Check whether the gate fits in each of the latest group of layers.
            # Start from the last of the already formed layers and go backwards, finding
            # the latest one where the pending gate doesn't fit: the _next_ one will be the
            # layer we need to put `g` in.

            idx = 0  # starting point for the iteration (the "best possible" scenario)
            @debug "Gate sites: " * join(_qbittag.(g_sites), " ")
            for (j, gl) in Iterators.reverse(enumerate(circ))
                actual_layer_idx = j
                @debug "Checking layer $actual_layer_idx. Free spots: " *
                    join(_qbittag.(freesites(circ, gl)), " ")
                if !issubset(g_sites, freesites(circ, gl))
                    @debug "No space available on this layer."
                    idx = j
                    break
                end
                # else the pending gate can "pass through" the current layer and we go on
                # to check the previous layer, to see if it fits there.
            end

            if idx + 1 > depth(circ)
                # This means that the first "blocking" gate is already the last one, so
                # we need to add a new layer to the circuit.
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
            eval(Meta.parse(new_gate_method))
        elseif line isa OpenQASM.Types.Barrier
            @warn "Barriers are not (yet) implemented. This line will be skipped."
        end
    end

    return circ
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
