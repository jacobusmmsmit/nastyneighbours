using StaticArrays
using Random
using CairoMakie
using JLD2
using ProgressMeter

using Coevolution

const strategies = SVector{16,SVector{4,Bool}}(
    SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
    SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
    SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
    SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
)

const group_agnostic_strategies = SVector{4,Int64}(findall(v -> (v[1] == v[3] && v[2] == v[4]), [SVector{4,Bool}((i - 1) >> shift & 1 != 0 for shift in 0:3) for i in 1:16]))

begin
    Z_group = 50
    Zs = (Z_group,)
    β = 1.0
    μ_s = 1 / 100 * sum(Zs)
    μ_g = 1 / 100 * sum(Zs)
    ξ = 0.1
    α = 0.9 # Assortment of interactions
    γ = 1.0 # Assortment of reproduction
    c = 1.0 # Cost of contribution
    ϵ_p = 0.01 # Error rate of production
    ϵ_c = 0.01 # Error rate of competition
    N = 100_000
end

S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
a = 1.0
pots = (SA[1+2., 2(1+2)], SA[0., 2(1+2)])
rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
@b main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)

l = 41
a_range = range(0, 4, length=l)
b_range = range(0, 4, length=l)

# In-group strategy
mean_strategy_count_matrix_grouped = [zeros(16) for i in b_range, j in a_range]
iterator = collect(Iterators.product(b_range, a_range))
@showprogress Threads.@threads for ij in 1:l^2
    b, a = iterator[ij]
    S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
    pots = (SA[1+b, 2(1+b)], SA[1+b, 2(1+b)])
    rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
    strategy_count_by_generation = main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)
    mean_strategy_count_matrix_grouped[ij] = dropdims(sum(strategy_count_by_generation[:, :, 9N÷10:end], dims=(1, 3)), dims=(1, 3)) ./ (N ÷ 10)
end

mean_strategy_count_matrix_agnostic = let
    M = [zeros(4) for i in b_range, j in a_range]
    for I in eachindex(mean_strategy_count_matrix_grouped)
        M[I] .= mean_strategy_count_matrix_grouped[I][group_agnostic_strategies]
    end
    M
end

# @save "data/mean_strategy_count_matrix_agnostic.jld2" mean_strategy_count_matrix_agnostic
# @load "data/mean_strategy_count_matrix_agnostic.jld2" mean_strategy_count_matrix_agnostic # v1 = b+1, v2 = 2(b+1)

begin
    fig = Figure(; size=(500, 450))
    for i in 1:2, j in 1:2
        fig[i, j] = GridLayout()
    end
    axs = [Axis(fig[i, j][1, 1], aspect=1, xticks=0:5, yticks=0:5, xlabel="Benefit (b)", ylabel = "Cost of Aggression (a)") for i in 1:2, j in 1:2]
    cmaps = map(color -> cgrad([:white, color]), strat_colours)
    iterator = ((2, 1), (1, 1), (2, 2), (1, 2))
    titles = ["Claim but don't produce", "Neither claim or produce", "Produce and claim", "Produce but don't claim"]
    for idx in 1:4
        i, j = iterator[idx]
        axs[idx].title = titles[idx]
        hm = heatmap!(
            axs[i, j], 
            b_range,
            a_range,
            colormap=cmaps[idx],
            getindex.(mean_strategy_count_matrix_agnostic, idx), 
            colorrange=(0, sum(Zs)),
        )
        Colorbar(fig[i, j][1, 2], hm)
    end
    for filetype in ["pdf", "png"]
        save("figures/revised/sim_2b+2.$filetype", fig)
    end
    fig
end

