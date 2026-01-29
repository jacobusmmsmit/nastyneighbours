struct RevisedParameters{N}
    Zs::NTuple{N,Int}
    β::Float64 # Strength of selection
    μ::Float64 # Mutation rate
    α::Float64 # Assortment of interactions
    γ::Float64 # Assortment of reproduction
    c::Float64 # Cost of contribution
    vs::SVector{4,Float64} # Vulnerability of each strategy
    as::SVector{4,Float64} # Cost of aggressing against each strategy
    pots::NTuple{2,SVector{3,Float64}}
    ϵ_p::Float64 # Error rate of production
    ϵ_c::Float64 # Error rate of competition
    n_migrants::Int64 # Migration rate
end

function average_utility(si::SVector{4,Bool}, gi::Integer, S, rp::RevisedParameters)
    (; Zs, α, c, vs, as, pots, ϵ_p, ϵ_c) = rp
    U_si = 0.0
    si_idx = evalpoly(2, si) + 1
    for (sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in, gj) in Iterators.product(false:true, false:true, false:true, false:true, 1:length(Zs))
        sj = SA[sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in]
        sj_idx = evalpoly(2, (sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in)) + 1
        sj_popsize = S[gj, sj_idx]
        sj_popsize == 0 && continue
        Z_i = Zs[gi]
        adjustment = (si == sj) & (gi == gj)
        if gi == gj
            sj_subset = SVector{2,Bool}(sj_claim_in, sj_produce_in)
            si_subset = SVector{2,Bool}(si[3], si[4])
            w = α
            pot = pots[1]
        else
            sj_subset = SVector{2,Bool}(sj_claim_out, sj_produce_out)
            si_subset = SVector{2,Bool}(si[1], si[2])
            w = (1 - α) / (length(Zs) - 1)
            pot = pots[2]
        end
        sj_weight = w * (sj_popsize - adjustment) / (α * Z_i + (1 - α) * (sum(Zs) - Z_i))
        if !(0 <= sj_weight <= 1)
            println("si: $si_idx, sj: $sj_idx, gi: $gi, gj: $gj, $adjustment, $sj_popsize")
            error("$si, $sj_popsize, sj_weight: $sj_weight, S_g: $(S[gi, :])")
        end
        payoff_matrix = get_payoff_matrix(pot, c, vs, as)
        payoff_matrix_with_errors = add_errors_to_payoff_matrix(payoff_matrix, ϵ_p, ϵ_c)
        U_si += sj_weight * payoff_from_interaction(si_subset, sj_subset, payoff_matrix_with_errors)
    end
    return U_si
end

function rand_S_initial_revised(Zs; strategy_set=1:16)
    res = zeros(Int, length(Zs), 16)
    rand_S_initial_revised!(res, Zs; strategy_set)
    return res
end

function rand_S_initial_revised!(res, Zs; strategy_set=1:16)
    res .= 0
    for (i, Z) in enumerate(Zs)
        for _ in 1:Z
            res[i, rand(strategy_set)] += 1
        end
    end
    return nothing
end

function main_simulation_loop(S_initial::AbstractMatrix{I}, N, rp::RevisedParameters{NG}; strategy_set=1:16) where {I<:Integer,NG}
    (; Zs, μ, β, n_migrants, γ) = rp
    n_agents = sum(Zs)
    S = copy(S_initial)
    T = zeros(Int, (NG, 16, N)) # Strategy by Group by Generation
    strategies = SVector{16,SVector{4,Bool}}(
        SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
        SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
        SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
        SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
    )
    group_weights = FrequencyWeights(dropdims(sum(S, dims=2), dims=2))
    migrants = zeros(Int, NG, n_migrants)
    swaps = Vector{Pair{Tuple{Int,Int},Tuple{Int,Int}}}(undef, n_migrants * NG ÷ 2)
    for G in 1:N
        for G_step in 1:n_agents # Strategy update
            # Sample two agents from the population to perform imitation
            gi = sample(1:NG, group_weights)
            if NG == 1
                gj = gi
            else
                gj = ifelse(rand() < γ, gi, rand_without(1:NG, gi))
            end
            i, j = sample_two_agents_without_replacement(S, gi, gj) # (i, j) .∈ Ref(1:16)
            if rand() < μ
                new_i = rand(strategy_set)
                S[gi, i] -= 1
                S[gi, new_i] += 1
            else
                # Calculate the utility of the strategies
                si = strategies[i]
                sj = strategies[j]
                U_i = average_utility(si, gi, S, rp)
                U_j = average_utility(sj, gj, S, rp)
                P_ij = inv(1 + exp(-β * (U_j - U_i)))
                # Imitate with probability P_ij which depends on utility difference.
                if rand() < P_ij
                    S[gi, i] -= 1
                    S[gi, j] += 1
                end
            end
        end
        # Migration: S[g, i] is the number of people playing the ith strategy in
        # the gth group. Sample n_migrants from each group.
        for g in 1:NG
            migrants[g, :] .= sample(1:Zs[g], n_migrants, replace=false)
        end
        # Pick two groups (potentially the same) with available migrants, choose
        # the first available migrants from those groups, add them to the list
        # of swaps.
        n_migrants_per_group = [n_migrants for _ in 1:NG]
        counters = [1 for _ in 1:NG]
        swap_idx = 0
        while sum(n_migrants_per_group) > 0
            swap_idx += 1
            # group_i and group_j may be the same
            group_i, group_j = sample_two_groups_without_replacement!(n_migrants_per_group)
            strat_i = get_strategy(migrants[group_i, counters[group_i]], @views(S[group_i, :]))
            strat_j = get_strategy(migrants[group_j, counters[group_j]], @views(S[group_j, :]))
            swaps[swap_idx] = Pair((group_i, strat_i), (group_j, strat_j))
            counters[group_i] += 1
            counters[group_j] += 1
        end
        # Perform the swaps
        for swap in swaps
            # [1] is the id, [2] is the strategy
            S[swap.first[1], swap.first[2]] -= 1
            S[swap.second[1], swap.second[2]] -= 1
            S[swap.first[1], swap.second[2]] += 1
            S[swap.second[1], swap.first[2]] += 1
        end
        @views T[:, :, G] .= S # Make a note of the current state of the population
    end
    return T
end

function rand_without(unitrange::UnitRange, j)
    _j = rand(unitrange.start:unitrange.stop-1)
    return ifelse(_j < j, _j, _j + 1)
end

function get_strategy(i, S_g)
    return findfirst(i .≤ cumsum(S_g))
end

function sample_two_groups_without_replacement!(n_migrants_per_group)
    cs = cumsum(n_migrants_per_group)
    i = rand(1:cs[end]) # choose a random agent
    group_i = findfirst(i .<= cs) # determine their group
    n_migrants_per_group[group_i] -= 1
    cumsum!(cs, n_migrants_per_group)
    j = rand(1:cs[end])
    group_j = findfirst(j .<= cs) # determine their group
    n_migrants_per_group[group_j] -= 1
    return group_i, group_j
end

function sample_two_agents_without_replacement(S, gi, gj) # 4 allocs
    S_gi = @views S[gi, :]
    if gi == gj
        cumtot = zero(MVector{16,Float64})
        tot = 0
        for idx in 1:16
            tot += S_gi[idx]
            cumtot[idx] = tot
        end
        i, j = sample(1:tot, 2, replace=true) # 2 allocations
        flag_i = false
        flag_j = false
        s_i = 0
        s_j = 0
        for idx in 1:16
            if !flag_i && i ≤ cumtot[idx]
                flag_i = true
                s_i = idx
            end
            if !flag_j && j ≤ cumtot[idx]
                flag_j = true
                s_j = idx
            end
            (flag_i && flag_j) && return (s_i, s_j)
        end
        error("Something went wrong: $i, $j, $tot, $cumtot")
    else
        S_gj = @views S[gj, :]
        si = sample(1:16, FrequencyWeights(S_gi))
        sj = sample(1:16, FrequencyWeights(S_gj))
        return (si, sj)
    end
end

"""
    get_cooperation_over_time(T)

Takes the output of a simulation and turns it into an `n_groups × 4 × n_timesteps`
array summarising how many individuals act cooperatively/competitively given
each of the four contexts.

Example output for one group, timestep pair: `[0, 5, 2, 1]` i.e. no one in this
group at this time cooperated with the outgroup, 5 cooperated with the ingroup,
2 competed with the outgroup, 1 competed with the ingroup.
"""
function get_cooperation_over_time(T)
    n_groups, n_strategies, n_timesteps = size(T)
    cooperation_per_context_per_group_per_timestep = zeros(eltype(T), n_groups, 8, n_timesteps)
    for timestep in 1:n_timesteps
        @views S = T[:, :, timestep]
        for (row_idx, strategy_counts) in enumerate(eachrow(S))
            for (strategy, strategy_count) in zip(strategies, strategy_counts)
                out_coop, out_comp, in_coop, in_comp = strategy
                out_strategy = evalpoly(2, SA[out_coop, out_comp]) + 1
                in_strategy = evalpoly(2, SA[in_coop, in_comp]) + 5
                cooperation_per_context_per_group_per_timestep[row_idx, out_strategy, timestep] += strategy_count
                cooperation_per_context_per_group_per_timestep[row_idx, in_strategy, timestep] += strategy_count
            end
        end
    end
    return cooperation_per_context_per_group_per_timestep
end