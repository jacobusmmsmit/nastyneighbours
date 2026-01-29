using StaticArrays
using Random
using CairoMakie

using Coevolution

begin
    Zs = (100, 10)
    β = 1.0
    μ = 1 / sum(Zs)
    α = 0.9 # Assortment of interactions
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
    n_repeats = 3
end

rp = RevisedParameters(Zs, β, μ, α, γ, c, vs, as, pots, ϵ_p, ϵ_c, n_migrants)
total_strategy_count_by_generation = sum(main_simulation_loop(rand_S_initial_revised(Zs), N, rp) for _ in 1:n_repeats)
mean_strategy_count_by_generation = total_strategy_count_by_generation ./ n_repeats

let
    fig = Figure(size=(600, 400))
    # Rows = [group 1, group 2,..., group N], columns = [out-group, in-group] 
    axs = [Axis(fig[by_group, to_relation], xlabel="Timestep", ylabel="Group $by_group", title="$(to_relation==1 ? "Out-group" : "In-group")") for (by_group, to_relation) in Iterators.product(1:length(Zs), 1:2)]
    cooperation_per_context_per_group_per_timestep = get_cooperation_over_time(mean_strategy_count_by_generation)
    for group in 1:length(Zs)
        engagement_per_context_per_timestep = @views cooperation_per_context_per_group_per_timestep[group, :, :]
        for (entry, is_to_ingroup, strategy) in zip(1:8, (1, 1, 1, 1, 2, 2, 2, 2), (1, 2, 3, 4, 1, 2, 3, 4))
            ax = axs[group, is_to_ingroup]
            cooperation_per_timestep = @views engagement_per_context_per_timestep[entry, :]
            lines!(ax, cooperation_per_timestep, color=Cycled(strategy))
        end
    end
    Legend(fig[1:2, 3], [LineElement(color=Cycled(i)) for i in 1:4], ["Free riders", "Producers", "Claimers", "PCs"])
    fig
end

Z2s_range = 10:10:100

T = main_simulation_loop(rand_S_initial_revised(Zs), N, rp)
get_strategic_diversity(T[:, :, :])
S = T[:, :, end]
mapslices(std, S ./ sum(S, dims=2), dims=2)

function get_strategic_diversity(T)
    # T is n_groups by 16 by n_generations
    strategic_diversity_over_time = mapslices(T, dims=(1, 2)) do S
        v = sum(S, dims=1)
        std(v / sum(v))
    end
    strategic_diversity_over_time = dropdims(strategic_diversity_over_time, dims=(1, 2))
end

function how_different_are_the_groups(T)
    v = mapslices(T, dims=(1, 2)) do S
        sum(abs, (S[1, :] ./ 100) .- (S[2, :] ./ sum(S[2, :]))) / 2
    end
    return dropdims(v, dims=(1, 2))
end


difference = map(Z2s_range) do (Z2)
    println("$Z2")
    Zs = (100, Z2)
    S_initial = rand_S_initial_revised(Zs)
    as = SA[a, a, a, a] # Cost of aggressing against each strategy
    pots = (SA[0, m_in*c, 2m_in*c], SA[0, m_out*c, 2m_out*c])
    rp = RevisedParameters(Zs, β, μ, α, γ, c, vs, as, pots, ϵ_p, ϵ_c, n_migrants)
    sum(1:n_repeats) do _
        T = main_simulation_loop(S_initial, N, rp)
        # mean(get_strategic_diversity(T))
        mean(how_different_are_the_groups(T))
    end / n_repeats
end

begin
    fig = Figure(size=(500, 300))
    ax = Axis(fig[1, 1], xlabel="Size of Group 2 (Group 1 has 100 Agents)", ylabel="Strategic Similarity between Groups")
    ylims!(ax, (0 - 0.05, 1 + 0.05))
    lines!(Z2s_range, 1 .- difference, linewidth=2, color=Cycled(2))    
    for filetype in ["png", "pdf"]
        save("./figures/revised/groupsize_nmigrants3.$filetype", fig)
    end
    display(fig)
end