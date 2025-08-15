using Random
using StaticArrays
using StatsBase

include("payoff_from_interaction.jl")

struct UnstructuredParameters
    Z::Int
    β::Float64
    μ::Float64
    c_p::Float64
    c_c::Float64
    m::Float64
    ϵ_p::Float64
    ϵ_c::Float64
end

function average_utility(si::SVector{2,Bool}, S, up::UnstructuredParameters)
    (; Z, c_p, c_c, m, ϵ_p, ϵ_c) = up
    U_si = 0.0
    for (sj_popsize, (sj_claim, sj_produce)) in zip(S, Iterators.product(false:true, false:true))
        # sj_popsize == 0 && continue
        sj = SA[sj_claim, sj_produce]
        adjustment = si == sj
        sj_weight = (sj_popsize - adjustment) / (Z - 1)
        U_si += sj_weight * payoff_from_interaction(si, sj, c_p, c_c, m, ϵ_p, ϵ_c)
    end
    return U_si
end

function rand_S_initial_unstructured(Z)
    res = zeros(Int, 4)
    rand_S_initial_unstructured!(res, Z)
    return res
end

function rand_S_initial_unstructured!(res, Z)
    res .= 0
    for _ in 1:Z
        res[rand(1:4)] += 1
    end
    return nothing
end

function main_simulation_loop(S_initial::AbstractVector{I}, N, up::UnstructuredParameters) where {I<:Integer}
    (; Z, μ, β) = up
    S = S_initial
    T = zeros(Int, (4, N)) # Strategy by Agents by Generation
    strategies = SVector{4,SVector{2,Bool}}(SA[0, 0], SA[1, 0], SA[0, 1], SA[1, 1])
    strategy_weights = FrequencyWeights(S, Z)
    for G in 1:N
        for _ in 1:Z
            i = sample(1:4, strategy_weights)
            si = strategies[i]
            if rand() < μ
                S[i] -= 1
                S[rand(1:4)] += 1
            else
                j = sample(1:4, strategy_weights)
                sj = strategies[j]
                U_i = average_utility(si, S, up)
                U_j = average_utility(sj, S, up)
                P_ij = inv(1 + exp(-β * (U_j - U_i)))
                if rand() < P_ij
                    S[i] -= 1
                    S[j] += 1
                end

            end
        end
        @views T[:, G] .= S
    end
    return T
end
