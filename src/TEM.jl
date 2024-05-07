module TEM

using ITensors, LindbladVectorizedTensors

export gate, GateName
include("gatename.jl")
#include("qelib1_gates_2.0.jl")
include("qelib1_gates_3.0.jl")
include("qelib1_gates_3.0vec.jl")

export qbit_registers, qbit_sites, qbit_map, interpret, gates
include("qasm_itensors_parser.jl")

include("spl_noise_model.jl")

end
