struct HierarchyParameters{N}
    Zs::NTuple{N,Int}
    β::Float64 # Strength of selection
    μ::Float64 # Mutation rate
    c_p::Float64
    c_c::Float64
    m_in::Float64
    m_out::Float64
    α::Float64 # Assortment
    ϵ_p::Float64 # Error rate of production
    ϵ_c::Float64 # Error rate of competition
    n_migrants::Int64 # Migration rate
    gini_coefficient::Float64
    hierarchy_strength::Float64
end

function HierarchyParameters(Zs::NTuple{N,Int}, β, μ, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c, n_migrants, gini, hierarchy_strength) where {N}
    return HierarchyParameters{N}(Zs, β, μ, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c, n_migrants, gini, hierarchy_strength)
end

function average_utility(si::SVector{4,Bool}, gi::Integer, S, hp::HierarchyParameters)
    (; Zs, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c) = hp
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
            m = m_in
        else
            sj_subset = SVector{2,Bool}(sj_claim_out, sj_produce_out)
            si_subset = SVector{2,Bool}(si[1], si[2])
            w = (1 - α) / (length(Zs) - 1)
            m = m_out
        end
        sj_weight = w * (sj_popsize - adjustment) / (α * Z_i + (1 - α) * (sum(Zs) - Z_i))
        if !(0 <= sj_weight <= 1)
            println("si: $si_idx, sj: $sj_idx, gi: $gi, gj: $gj, $adjustment, $sj_popsize")
            error("$si, $sj_popsize, sj_weight: $sj_weight, S_g: $(S[gi, :])")
        end
        U_si += sj_weight * payoff_from_interaction(si_subset, sj_subset, c_p, c_c, m, ϵ_p, ϵ_c)
    end
    return U_si
end

function calculate_hierarchy_bonus(average_utilities, S_g, hp)
    avg_utility_in_group = sum(S_g[i] * average_utilities[i] for i in 1:16) / sum(S_g)
    residuals = average_utilities .- avg_utility_in_group
    return SVector{16,Float64}(ifelse(r > 0, r^hp.inequality, -r^hp.inequality) for r in residuals)
    # return SVector{16,Float64}(average_utilities[i]^hp.inequality / sum(S_g[i] * average_utilities[i]^hp.inequality) for i in 1:16)
end

function rand_S_initial_hierarchy(Zs)
    res = zeros(Int, length(Zs), 16)
    rand_S_initial_hierarchy!(res, Zs)
    return res
end

function rand_S_initial_hierarchy!(res, Zs)
    res .= 0
    for (i, Z) in enumerate(Zs)
        for _ in 1:Z
            res[i, rand(1:16)] += 1
        end
    end
    return nothing
end

function main_simulation_loop(S_initial::AbstractMatrix{I}, N, hp::HierarchyParameters{NG}) where {I<:Integer,NG}
    (; Zs, μ, β, n_migrants) = hp
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
    swaps = Vector{Pair{Tuple{Int, Int}, Tuple{Int, Int}}}(undef, n_migrants*NG÷2)
    for G in 1:N
        for G_step in 1:n_agents # Strategy update
            # Sample two agents in the same group as imitation only happens within the group
            gi = sample(1:NG, group_weights)
            S_gi = @views S[gi, :]
            i, j = sample_two_agents_without_replacement(S_gi) # (i, j) .∈ Ref(1:16)
            # Calculate the utility of all in-group members
            average_utilities = SVector{16,Float64}(average_utility(s, gi, S, hp) for s in strategies)
            redistributed_utilities = redistribute_utilities(average_utilities, S_gi, hp)
            # println("average utilities: $(round.(average_utilities; sigdigits=2))")
            # println("redistributed uts: $(round.(redistributed_utilities; sigdigits=2))")
            # println("         n agents: $(round.(S_gi; sigdigits=2))")
            if rand() < μ
                new_i = rand(1:16)
                S[gi, i] -= 1
                S[gi, new_i] += 1
            else
                U_i = redistributed_utilities[i]
                U_j = redistributed_utilities[j]
                P_ij = inv(1 + exp(-β * (U_j - U_i)))
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
        # Pick two groups with available migrants, choose the first available
        # migrants from those groups, add them to the list of swaps.
        n_migrants_per_group = [n_migrants for _ in 1:NG]
        counters = [1 for _ in 1:NG]
        swap_idx = 0
        while sum(n_migrants_per_group) > 0
            swap_idx += 1
            group_i, group_j = sample_two_groups_without_replacement!(n_migrants_per_group)
            strat_i = get_strategy(migrants[group_i, counters[group_i]], @views(S[group_i, :]))
            strat_j = get_strategy(migrants[group_j, counters[group_j]], @views(S[group_j, :]))
            swaps[swap_idx] = Pair((group_i, strat_i), (group_j, strat_j))
            counters[group_i] += 1
            counters[group_j] += 1
        end
        # Perform the swaps
        for swap in swaps
            S[swap.first[1], swap.first[2]] -= 1
            S[swap.second[1], swap.second[2]] -= 1
            S[swap.first[1], swap.second[2]] += 1
            S[swap.second[1], swap.first[2]] += 1
        end
        @views T[:, :, G] .= S # Make a note of the current state of the population
    end
    return T
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


function redistribute_utilities(average_utilities, S_gi, hp)
    au_nz = average_utilities[S_gi .> 0]
    S_nz = S_gi[S_gi .> 0]
    v = sortperm(au_nz)
    iv = invperm(v)
    n = S_nz[v]
    y = au_nz[v]
    y_new, _, _, _, _ = interpolate_lorenz_bins(n, y, hp.gini_coefficient, hp.hierarchy_strength)
    # println("avg: $(round.(average_utilities, sigdigits=2))")
    # println("$(round.(y_new[iv], sigdigits=2))")
    out = zeros(16)
    out[S_gi .> 0] .= y_new[iv]
    out[S_gi .== 0] .= average_utilities[S_gi .== 0]
    return out
end

function sample_two_agents_without_replacement(S_g) # 4 allocs
    cumtot = zero(MVector{16,Float64})
    tot = 0
    for idx in 1:16
        tot += S_g[idx]
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
end

# Some AI generated helper functions
# --- 1. Lorenz interpolation function ---
function interpolate_lorenz_bins(n::AbstractVector, y::AbstractVector, G_target::Real, λ::Real)
    # --- Compute empirical Lorenz curve ---
    N = sum(n)
    Y = sum(n .* y)
    cum_n = cumsum(n)
    cum_inc = cumsum(n .* y)
    
    p_emp = vcat(0.0, cum_n ./ N)
    L_emp = vcat(0.0, cum_inc ./ Y)
    
    # --- Synthetic power-Lorenz curve ---
    if !(0 <= G_target < 1)
        error("G_target must be in [0,1).")
    end
    k = (1 + G_target) / (1 - G_target)
    L_synth = p_emp .^ k
    
    # --- Interpolate ---
    L_mix = (1 - λ) .* L_emp .+ λ .* L_synth

    L_mix[1] ≈ 0.0 || error("$L_mix, $n")
    # L_mix[end] ≈ 1.0 || error("$L_mix, $n")
    # L_mix = setindex(L_mix, 1, length(n))
    L_mix[end] = 1.0
    
    # --- Compute constant per-bin incomes ---
    y_new = zeros(length(n))
    for i in 1:length(n)
        L_start = interp_linear(p_emp[i], p_emp, L_mix)
        L_end   = interp_linear(p_emp[i+1], p_emp, L_mix)
        bin_share = L_end - L_start
        y_new[i] = bin_share * Y / n[i]
    end
    
    return y_new, p_emp, L_emp, L_synth, L_mix
end

# --- Helper: linear interpolation ---
function interp_linear(x, xp::AbstractVector, yp::AbstractVector)
    if x <= xp[1]
        return yp[1]
    elseif x >= xp[end]
        return yp[end]
    else
        idx = searchsortedfirst(xp, x)
        x0, x1 = xp[idx-1], xp[idx]
        y0, y1 = yp[idx-1], yp[idx]
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    end
end