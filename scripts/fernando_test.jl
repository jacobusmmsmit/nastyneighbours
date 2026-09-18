using StaticArrays
using Random
using CairoMakie
using JLD2
using ProgressMeter
using Interpolations
using LinearAlgebra: norm
using UniformStreamlines

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
    μ_s = 1 / (1 * sum(Zs))
    μ_g = 1 / (1 * sum(Zs))
    ξ = 1.0 # Likelihood of strategic update as opposed to a group update
    α = 1.0 # Assortment of interactions
    γ = 1.0 # Assortment of reproduction
    c = 1.0 # Cost of contribution
    ϵ_p = 0.01 # Error rate of production
    ϵ_c = 0.01 # Error rate of competition
    N = 10 * Z_group^2
end

v1(b) = 1 + b
v2(b) = 2(1 + b)
b = 1.0
a = 1.0
pots = (SA[v1(b), v2(b)], SA[v1(b), v2(b)])

S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
@b main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)

begin
    Profile.Allocs.clear()
    main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)
    Profile.Allocs.@profile sample_rate=1 main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)
    PProf.Allocs.pprof(from_c=false)
end

l = 6
α_range = range(0, 1, length=l)
γ_range = range(0, 1, length=l)

# Collect only the four agnostic strategies
mean_strategy_count_matrix_grouped = let
    M = [zeros(4) for i in α_range, j in γ_range]
    iterator = collect(Iterators.product(α_range, γ_range))
    @showprogress Threads.@threads for ij in 1:(l^2)
        α, γ = iterator[ij]
        S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
        rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
        strategy_count_by_generation = main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)
        burn_in_period = N ÷ 10
        collection_period = N - burn_in_period
        M[ij] = dropdims(sum(strategy_count_by_generation[:, :, (end-collection_period+1):end], dims=(1, 3)), dims=(1, 3))[[1, 6, 11, 16]] ./ collection_period
    end
    M
end

begin
    fig2 = Figure(; size=(600, 600))
    ga = fig2[1, 1] = GridLayout()
    # gb = fig2[1, 2] = GridLayout()
    begin
        for i in 1:2, j in 1:2
            idx = 2(i - 1) + (j - 1)
            ax = Axis(ga[3-j, 2i-1]; aspect=1)
            ax.xlabel = "Interaction Assortment (α)"
            ax.ylabel = "Reproduction Assortment (γ)"
            ax.title = labels[idx+1]
            hm = heatmap!(
                ax,
                α_range,
                γ_range,
                getindex.(mean_strategy_count_matrix_grouped, idx + 1),
                colorrange=(0, sum(Zs)),
                colormap=cgrads[idx+1],
            )
            cb = Colorbar(ga[3-j, 2i], hm, label="Number of agents", tellheight=true)
            cb.height = Relative(0.73)
            # limits!(ax, (0, 4), (0, 4))
        end
        for (label, pos) in zip(["a", "b", "c", "d"], [[1, 1], [1, 3], [2, 1], [2, 3]])
            Label(ga[pos[1], pos[2], TopLeft()], label,
                fontsize=26,
                font=:bold,
                padding=(0, 5, 0, 0),
                halign=:right)
        end

    end
    display(fig2)
end