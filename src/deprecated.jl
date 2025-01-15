"""
    gatelayers(gates::Vector{ITensor})

Return a list of MPOs, each one representing a layer of the circuit.
A layer is created multiplying together as many adjacent gates as possible, i.e. the
construction of the MPO ends as soon as the next gate in line acts on a qbit which is
already acted on by the gates already in the MPO.
This way, all gates in a layer commute with each other and can be executed in any order.
"""
function gatelayers(gates::Vector{ITensor})
    Base.depwarn(
        "The `gatelayers(gates)` function will be deprecated in a future release. " *
        "To organize the gates into layers, please use the " *
        "`parsecircuit(code)` instead, loading the circuit directly from an OpenQASM " *
        "program into a QuantumCircuit object. The new " *
        "also uses a \"greedier\" " *
        "layering algorithm that should produce fewer layers.",
        :gatelayers,
    )
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
