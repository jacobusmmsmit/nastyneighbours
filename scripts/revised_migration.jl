using StaticArrays
using Random
using CairoMakie

using Coevolution

begin
    Zs = (100, 100)
    β = 1.0
    μ = 1 / sum(Zs)
    α = 0.5 # Assortment of interactions
    γ = 1.0 # Assortment of reproduction
    vs = SA[1.0, 1, 1, 1] # Vulnerability of each strategy
    c = 1.0 # Cost of contribution
    m_in = 2.5
    m_out = 2.0
    a = 1.5
    as = SA[a, a, a, a] # Cost of aggressing against each strategy
    pots = (SA[0, m_in*c, 2m_in*c], SA[0, m_out*c, 2m_out*c])
    ϵ_p = 0.01 # Error rate of production
    ϵ_c = 0.01 # Error rate of competition
    n_migrants = 10 # Migration rate
    N = 1000
    n_repeats = 10
end

S_initial = rand_S_initial_revised(Zs)
rp = RevisedParameters(Zs, β, μ, α, γ, c, vs, as, pots, ϵ_p, ϵ_c, n_migrants)
total_strategy_count_by_generation = sum(main_simulation_loop(rand_S_initial_revised(Zs), N, rp) for _ in 1:n_repeats)
mean_strategy_count_by_generation = total_strategy_count_by_generation ./ n_repeats


let
    fig = Figure(size=(600, 100 + 150 * length(Zs)))
    # Rows = [group 1, group 2,..., group N], columns = [out-group, in-group] 
    axs = [Axis(fig[by_group, to_relation], xlabel="Timestep", ylabel="Group $by_group", title="$(to_relation==1 ? "Out-group" : "In-group")") for (by_group, to_relation) in Iterators.product(1:length(Zs), 1:2)]
    cooperation_per_context_per_group_per_timestep = get_cooperation_over_time(mean_strategy_count_by_generation)
    for group in 1:length(Zs)
        engagement_per_context_per_timestep = @views cooperation_per_context_per_group_per_timestep[group, :, :]
        display(engagement_per_context_per_timestep)
        for (entry, is_to_ingroup, strategy) in zip(1:8, (1, 1, 1, 1, 2, 2, 2, 2), (1, 2, 3, 4, 1, 2, 3, 4))
            ax = axs[group, is_to_ingroup]
            cooperation_per_timestep = @views engagement_per_context_per_timestep[entry, :]
            lines!(ax, cooperation_per_timestep, color=Cycled(strategy))
        end
    end
    linkaxes!(axs...)
    Legend(fig[1:length(Zs), 3], [LineElement(color=Cycled(i)) for i in 1:4], ["Free riders", "Producers", "Claimers", "PCs"])
    fig
end

n_migrants_range = 0:1:minimum(Zs)

S_initial = rand_S_initial_revised(Zs)
difference = map(n_migrants_range) do (n_migrants)
    println("$n_migrants")
    rand_S_initial_revised!(S_initial, Zs)
    as = SA[a, a, a, a] # Cost of aggressing against each strategy
    pots = (SA[0, m_in*c, 2m_in*c], SA[0, m_out*c, 2m_out*c])
    rp = RevisedParameters(Zs, β, μ, α, γ, c, vs, as, pots, ϵ_p, ϵ_c, n_migrants)
    sum(1:n_repeats) do _
        T = main_simulation_loop(S_initial, N, rp)
        mapsl = mapslices(msc -> sum(abs, msc[1, :] .- msc[2, :]) ÷ 2, T, dims=(1, 2))
        res = mean(dropdims(mapsl, dims=(1, 2)))
    end / n_repeats
end

begin
    fig = Figure(size=(500, 300))
    ax = Axis(fig[1, 1], xlabel="Number of Migrants per Generation", ylabel="Strategic Similarity between Groups")
    # ylims!(ax, (0-0.05, 1+0.05))
    lines!(n_migrants_range, (Zs[1] .- difference) / Zs[1], linewidth=2)
    display(fig)
    save("./figures/revised/migration.pdf", fig)
end