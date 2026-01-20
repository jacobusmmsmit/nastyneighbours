struct StructuredParameters
    Z_1::Int
    Z_2::Int
    β::Float64
    μ::Float64
    c_p::Float64
    c_c::Float64
    m_in::Float64
    m_out::Float64
    ϵ_p::Float64
    ϵ_c::Float64
    α::Float64
end

struct StructuredParameters{F}
    Z_1::Int
    Z_2::Int
    β::Float64
    μ::Float64
    c::Float64
    vs::SVector{4,Float64}
    as::SVector{4,Float64}
    pot_function::F #F: pot parameters (b0, b1, b2) → SVector{3, Float64}
    ϵ_p::Float64
    ϵ_c::Float64
    α::Float64
end

function average_utility(si::SVector{4,Bool}, gi::Bool, S, sp::StructuredParameters)
    (; Z_1, Z_2, c, vs, as, m_in, m_out, ϵ_p, ϵ_c, α) = sp
    U_si = 0.0 # TODO: Test this function
    for (sj_popsize, (sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in, gj)) in zip(S, Iterators.product(false:true, false:true, false:true, false:true, false:true))
        # sj_popsize == 0 && continue
        Z_i = (Z_1, Z_2)[gi+1]
        Z_noti = (Z_2, Z_1)[gi+1]
        sj = SA[sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in]
        adjustment = (si == sj) & (gi == gj)
        if gi == gj
            sj_subset = SVector{2,Bool}(sj_claim_in, sj_produce_in)
            si_subset = SVector{2,Bool}(si[3], si[4])
            w = α
            m = m_in
        else
            sj_subset = SVector{2,Bool}(sj_claim_out, sj_produce_out)
            si_subset = SVector{2,Bool}(si[1], si[2])
            w = 1 - α
            m = m_out
        end
        sj_weight = w * (sj_popsize - adjustment) / (α * (Z_i - 1) + (1 - α) * Z_noti)
        if !(0 <= sj_weight <= 1)
            println(S[1:16])
            println(sj_popsize)
            error("$si, $sj_popsize, sj_weight: $sj_weight") # = (1 - $α) * ($sj_popsize - $adjustment) / ($α*($Z_i - 1) + (1-$α)*$Z_noti)")
        end
        pot::SVector{3,Float64} = sp.pot_function(b0, m * c, m * c)
        payoff_matrix = get_payoff_matrix(pot, c, vs, as)
        pmwe = add_errors_to_payoff_matrix(payoff_matrix, ϵ_p, ϵ_c)
        U_si += sj_weight * payoff_from_interaction(si_subset, sj_subset, pmwe)
    end
    return U_si
end

function rand_S_initial_structured(Z_1, Z_2)
    res = zeros(Int, 32)
    rand_S_initial_structured!(res, Z_1, Z_2)
    return res
end

function rand_S_initial_structured!(res, Z_1, Z_2)
    res .= 0
    for _ in 1:Z_1
        res[rand(1:16)] += 1
    end
    for _ in 1:Z_2
        res[rand(17:32)] += 1
    end
    return nothing
end

function main_simulation_loop(S_initial::AbstractVector{I}, N, sp::StructuredParameters) where {I<:Integer}
    (; Z_1, Z_2, μ, β) = sp
    S = copy(S_initial)
    T = zeros(Int, (32, N)) # Strategy by Agents by Generation
    strategies = SVector{32,SVector{4,Bool}}(
        SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
        SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
        SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
        SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
        SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
        SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
        SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
        SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1]
    )
    strategy_weights = FrequencyWeights(S, Z_1 + Z_2)
    for G in 1:N
        for _ in 1:(Z_1+Z_2)
            i = sample(1:32, strategy_weights)
            gi = i > 16
            si = strategies[i]
            if rand() < μ
                S[i] -= 1
                new_i = rand(1:16) + 16gi
                abs(i - new_i) > 15 && error("$i and $new_i not in same group")
                S[new_i] += 1
            else
                # Imitation only happens within group
                in_group_weights = @views strategy_weights[(1:16).+16gi]
                j = sample((1:16) .+ 16gi, FrequencyWeights(in_group_weights, (Z_1, Z_2)[gi+1]))
                sj = strategies[j]
                U_i = average_utility(si, gi, S, sp)
                U_j = average_utility(sj, gi, S, sp)
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
