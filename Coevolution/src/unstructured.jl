struct UnstructuredParameters
    Z::Int
    β::Float64
    μ::Float64
    c::Float64
    vs::SVector{4, Float64}
    as::SVector{4, Float64}
    pot::SVector{3, Float64}
    ϵ_p::Float64
    ϵ_c::Float64
end

function average_utility(si::SVector{2, Bool}, S, payoff_matrix, Z)
    U_si = 0.0
    for (sj_popsize, (sj_claim, sj_produce)) in zip(S, Iterators.product(false:true, false:true))
        # sj_popsize == 0 && continue
        sj = SA[sj_claim, sj_produce]
        adjustment = si == sj
        sj_weight = (sj_popsize - adjustment) / (Z - 1)
        U_si += sj_weight * payoff_from_interaction(si, sj, payoff_matrix)
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

function main_simulation_loop(S_initial::AbstractVector{I}, N, up::UnstructuredParameters) where {I <: Integer}
    (; Z, μ, β) = up
    payoff_matrix = get_payoff_matrix(up.pot, up.c, up.vs, up.as)
    payoff_matrix_with_errors = add_errors_to_payoff_matrix(payoff_matrix, up.ϵ_p, up.ϵ_c)
    S = S_initial
    T = zeros(Int, (4, N)) # Strategy by Agents by Generation
    strategies = SVector{4, SVector{2, Bool}}(SA[0, 0], SA[1, 0], SA[0, 1], SA[1, 1])
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
                U_i = average_utility(si, S, payoff_matrix_with_errors, Z)
                U_j = average_utility(sj, S, payoff_matrix_with_errors, Z)
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
