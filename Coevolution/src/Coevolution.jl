module Coevolution

using Random

using StaticArrays
using StatsBase
using CairoMakie


# Outgroup strategy is first two bits of strategy, in-group is third and fourth.
# "Produce" bits are 2 and 4, "Claim" bits are 1 and 3.
strategies = SVector{16,SVector{4,Bool}}(
    SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
    SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
    SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
    SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
)

export one_at, to_bin, repeat_vector
include("helper_functions.jl")

export strategic_errors, get_payoff_matrix, payoff_from_interaction, get_pot
include("payoff_from_interaction.jl")

export labels, strat_colours, cgrads
include("constants.jl")

export main_simulation_loop, average_utility # Common to both unstructured and structured
export UnstructuredParameters, rand_S_initial_unstructured, rand_S_initial_unstructured!
include("unstructured.jl")

export StructuredParameters, rand_S_initial_structured, rand_S_initial_structured!
include("structured.jl")

export RevisedParameters, rand_S_initial_revised, rand_S_initial_revised!, sample_two_agents_without_replacement, get_cooperation_over_time
include("revised.jl")

end # module Coevolution
