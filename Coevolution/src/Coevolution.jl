module Coevolution

using Random

using StaticArrays
using StatsBase
using CairoMakie

export one_at, to_bin, repeat_vector
include("helper_functions.jl")

export strategic_errors, get_payoff_matrix, payoff_from_interaction
include("payoff_from_interaction.jl")

export labels, strat_colours, cgrads
include("constants.jl")

export main_simulation_loop, average_utility # Common to both unstructured and structured
export UnstructuredParameters, rand_S_initial_unstructured, rand_S_initial_unstructured!
include("unstructured.jl")

export StructuredParameters, rand_S_initial_structured, rand_S_initial_structured!
include("structured.jl")

end # module Coevolution
