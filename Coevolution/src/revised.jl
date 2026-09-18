struct RevisedParameters{N}
    initial_Zs::NTuple{N,Int}
    β::Float64 # Strength of selection of strategies and groups
    μ_s::Float64 # Mutation rate of strategies
    μ_g::Float64 # Mutation rate of group tags
    ξ::Float64 # Probability of strategy update (versus a group tag update)
    α::Float64 # Assortment of interactions
    γ::Float64 # Assortment of reproduction
    c::Float64 # Cost of contribution
    a::Float64 # Cost of aggressing
    pots::NTuple{2,SVector{2,Float64}} # [v0, v1, v2] total utility from n agents cooperating
    ϵ_p::Float64 # Error rate of production
    ϵ_c::Float64 # Error rate of competition
end

function average_utility(si::SVector{4,Bool}, gi::Integer, S, rp::RevisedParameters{N}) where {N}
    (; α, c, a, pots, ϵ_p, ϵ_c) = rp
    U_si = 0.0
    si_idx = evalpoly(2, si) + 1
    # Loop over all possible partners by going over each strategy in each group.
    for (sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in, gj) in Iterators.product(false:true, false:true, false:true, false:true, 1:N)
        sj = SA[sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in] # Strategy of partner `j`
        # Convert the strategy `sj` from binary to decimal (add one for one-based indexing)
        sj_idx = evalpoly(2, (sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in)) + 1
        Z_j = S[gj, sj_idx] # Number of agents in group `gj` playing strategy `sj`
        # If agent `i` and `j` are in the same group, playing the same strategy,
        # then we must be careful that `i` and `j` may be referring to the same
        # individual. Thus we have to explicitly remove `i` from the pool of
        # partners by subtracting an `adjustment` (either 0 or 1).
        adjustment = (si == sj) && (gi == gj)
        # If no one (else) is playing this strategy (sj) in this group (gj) we can skip this iteration.
        Z_j - adjustment == 0 && continue
        Z_j - adjustment < 0 && error("What? Z_j = $Z_j, adjustment = $adjustment, but S[si, gi] = $(S[si_idx, gi]), dumping inputs: $si\n$gi\n$S\n$rp")
        # Count the number of people inside and outside of group gi, adjust for
        # not counting `i` itself.
        Z_i, Z_noti = let
            inside_total = 0
            for i in axes(S, 2)
                @inbounds inside_total += S[gi, i]
            end
            (inside_total - 1, sum(S) - inside_total - 1)
        end
        # Given that different bits in the strategy of each agent refer to their
        # actions towards different agents, either in-group or out-group, we
        # construct their "strategy subset" based on which bits are relevant.
        if gi == gj
            sj_subset = SVector{2,Bool}(sj_claim_in, sj_produce_in)
            si_subset = SVector{2,Bool}(si[3], si[4])
            w = α
            aggression_ratio = 1 // 2
            pot = pots[1]
        else
            sj_subset = SVector{2,Bool}(sj_claim_out, sj_produce_out)
            si_subset = SVector{2,Bool}(si[1], si[2])
            w = (1 - α) / (N - 1)
            aggression_ratio = let
                num = 0
                denom = 0
                for s in 2:2:16
                    num += S[gi, s]
                    denom += S[gi, s] + S[gj, s]
                end
                denom > 0 ? num // denom : 1 // 2
            end
            pot = pots[2]
        end
        sj_weight = w * (Z_j - adjustment) / (α * Z_i + (1 - α) * Z_noti)
        payoff_matrix = get_payoff_matrix(pot, c, a, aggression_ratio)
        payoff_matrix_with_errors = add_errors_to_payoff_matrix(payoff_matrix, ϵ_p, ϵ_c)
        U_si += sj_weight * payoff_from_interaction(si_subset, sj_subset, payoff_matrix_with_errors)
    end
    return U_si
end

function rand_S_initial_revised(initial_Zs; strategy_set=1:16)
    res = zeros(Int, length(initial_Zs), 16)
    rand_S_initial_revised!(res, initial_Zs; strategy_set)
    return res
end

function rand_S_initial_revised!(res, initial_Zs; strategy_set=1:16)
    res .= 0
    for (i, Z) in enumerate(initial_Zs)
        for _ in 1:Z
            res[i, rand(strategy_set)] += 1
        end
    end
    return nothing
end

function main_simulation_loop(S_initial::AbstractMatrix{I}, N, rp::RevisedParameters{NG}; strategy_set=1:16) where {I<:Integer,NG}
    # Setup phase
    (; initial_Zs, ξ, μ_s, μ_g, β, γ) = rp
    n_agents = sum(initial_Zs)
    S = copy(S_initial)
    T = zeros(Int, (NG, 16, N)) # Strategy by Group by Generation
    strategies = SVector{16,SVector{4,Bool}}(
        SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
        SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
        SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
        SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
    )
    # group_weight_vector = MVector{NG,I}(dropdims(sum(S, dims=2), dims=2))
    group_weight_vector = dropdims(sum(S, dims=2), dims=2)
    # Simulation loop phase
    for G in 1:N
        for G_step in 1:n_agents
            ## Strategy update or group update?
            # Sample two agents from the population to perform imitation
            gi = sample_group(group_weight_vector, n_agents)
            imitate_ingroup = nothing
            if NG == 1
                gj = gi
            else
                imitate_ingroup = rand() < γ
                gj = ifelse(
                    group_weight_vector[gi] == n_agents || (imitate_ingroup && (group_weight_vector[gi] > 1)),
                    gi,
                    sample_without(group_weight_vector, gi))
            end
            i, j = sample_two_agents_without_replacement(S, gi, gj) # (i, j) .∈ Ref(1:16)
            if rand() < ξ
                # Strategy update
                if rand() < μ_s
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
            else
                ## Group update
                if rand() < μ_g
                    # Randomly join a new group, this group may have members or not
                    new_g = rand(1:NG)
                    S[gi, i] -= 1
                    S[new_g, i] += 1
                    group_weight_vector[gi] -= 1
                    group_weight_vector[new_g] += 1
                else
                    # Join the group of another agent depending on their payoff.
                    # First, calculate the utility of the strategies.
                    si = strategies[i]
                    sj = strategies[j]
                    U_i = average_utility(si, gi, S, rp)
                    U_j = average_utility(sj, gj, S, rp)
                    P_ij = inv(1 + exp(-β * (U_j - U_i)))
                    # Then imitate group i.e. join that group with probability P_ij.
                    # Possible (depending on code) also imitate strategy.
                    if rand() < P_ij
                        S[gi, i] -= 1
                        S[gj, i] += 1 # j for mutate strategy, i for mutate group only
                        group_weight_vector[gi] -= 1
                        group_weight_vector[gj] += 1
                    end
                end
            end
        end
        @views T[:, :, G] .= S # Make a note of the current state of the population
    end
    return T
end

function rand_without(unitrange::UnitRange, j)
    _j = rand(unitrange.start:(unitrange.stop-1))
    return ifelse(_j < j, _j, _j + 1)
end

function get_strategy(i, S_g)
    return findfirst(i .≤ cumsum(S_g))
end

@inline function sample_without(weights, j, total)
    total_without = total - weights[j]
    r = rand(1:total_without)
    @inbounds for g in eachindex(weights)
        g == j && continue
        r -= weights[g]
        r ≤ 0 && return g
    end
    return lastindex(weights)
end

@inline function sample_two_same_group(S_gi)
    tot = sum(S_gi)

    i = rand(1:tot)
    j = rand(1:tot - 1)
    j += (j >= i)

    s_i = 1
    @inbounds while i > S_gi[s_i]
        i -= S_gi[s_i]
        s_i += 1
    end

    s_j = 1
    @inbounds while j > S_gi[s_j]
        j -= S_gi[s_j]
        s_j += 1
    end

    return s_i, s_j
end

@inline function sample_strategy(S, g)
    tot = 0
    @inbounds for idx in 1:16
        tot += S[g, idx]
    end
    r = rand(1:tot)
    @inbounds for idx in 1:16
        r -= S[g, idx]
        r ≤ 0 && return idx
    end
    return 16
end

@inline function sample_strategy(S, g, total)
    r = rand(1:total)

    @inbounds for i in 1:16
        r -= S[g, i]
        r ≤ 0 && return i
    end

    return 16
end

@inline function sample_two_same_group(S, g, total)
    i = rand(1:total)

    j = rand(1:total - 1)
    j += (j >= i)

    s_i = 1
    @inbounds while i > S[g, s_i]
        i -= S[g, s_i]
        s_i += 1
    end

    s_j = 1
    @inbounds while j > S[g, s_j]
        j -= S[g, s_j]
        s_j += 1
    end

    return s_i, s_j
end

function sample_two_agents_without_replacement(S, gi, gj) # 4 allocs
    S_gi = @views S[gi, :]
    if gi == gj
        si, sj = sample_two_same_group(S_gi)
    else
        si = sample_strategy(S, gi, tot_i)
        sj = sample_strategy(S, gj, tot_j)
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

@inline function sample_group(weights, total)
    r = rand(1:total)
    @inbounds for g in eachindex(weights)
        r -= weights[g]
        r ≤ 0 && return g
    end
    return lastindex(weights)
end