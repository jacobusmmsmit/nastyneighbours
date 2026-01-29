using StaticArrays
using Random
using CairoMakie
using StatsBase

using Coevolution

begin
    Zs = (30, 30, 30, 30)
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
    N = 200
    n_repeats = 50
end

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

ngs_range = (2, 4, 6, 8, 10, 12)
lcm_ngs = lcm(ngs_range...)

# TODO: Define a function that has equal "strategic diversity" for a given
# random initialisation for any number of groups.
function get_strategic_diversity(T)
    # T is n_groups by 16 by n_generations
    strategic_diversity_over_time = mapslices(T, dims=(1, 2)) do S
        v = sum(S, dims=1)
        std(v/sum(v))
    end
    strategic_diversity_over_time = dropdims(strategic_diversity_over_time, dims=(1, 2))
end

function get_average_strategic_diversity_within_group(T)
    # T is n_groups by 16 by n_generations
    _out = mapslices(T, dims=(1, 2)) do S
        v = mapslices(std, S./sum(S, dims=2), dims=2)
        mean(v)
    end
    return dropdims(_out, dims=(1, 2))
end


T = main_simulation_loop(rand_S_initial_revised(Zs), N, rp)
get_strategic_diversity(T)
get_average_strategic_diversity_within_group(T)

data = map(ngs_range) do NG
    Zs = Tuple(lcm_ngs ÷ NG for _ in 1:NG)
    S_initial = rand_S_initial_revised(Zs)
    as = SA[a, a, a, a] # Cost of aggressing against each strategy
    pots = (SA[0, m_in*c, 2m_in*c], SA[0, m_out*c, 2m_out*c])
    rp = RevisedParameters(Zs, β, μ, α, γ, c, vs, as, pots, ϵ_p, ϵ_c, n_migrants)
    mean_strategic_diversity = zeros(200)
    mean_strategic_diversity_in_group = zeros(200)
    for _ in 1:n_repeats
        T = main_simulation_loop(S_initial, N, rp)
        mean_strategic_diversity += get_strategic_diversity(T)
        mean_strategic_diversity_in_group += get_average_strategic_diversity_within_group(T)
    end
    mean_strategic_diversity ./ n_repeats, mean_strategic_diversity_in_group ./ n_repeats
end

begin
    fig = Figure(size=(500, 300))
    ax = Axis(fig[1, 1], xlabel = "Generation", ylabel = "Strategic Diversity of Population")
    for (v_idx, (v, _)) in enumerate(data)
        lines!(1:N, v, label="$(ngs_range[v_idx])", color = ngs_range[v_idx]/12, colorrange = (0, 1), linewidth=3)
    end
    Legend(fig[1, 2], [LineElement(color=ng/12, linewidth=3, colorrange=(0, 1)) for ng in reverse(ngs_range)], ["$ng" for ng in reverse(ngs_range)], "Number of Groups")
    # axislegend(ax, "N Groups", nbanks=3, position=:lt)
    for filetype in ["png", "pdf"]
        save("./figures/revised/ngroups_population_level_10migrants.$filetype", fig)
    end
    fig
end

begin
    fig = Figure(size=(500, 300))
    ax = Axis(fig[1, 1], xlabel = "Generation", ylabel = "Strategic Diversity within Groups")
    for (v_idx, (_, v)) in enumerate(data)
        lines!(1:N, v, label="$(ngs_range[v_idx])", color = ngs_range[v_idx]/12, colorrange = (0, 1), linewidth=3)
    end
    Legend(fig[1, 2], [LineElement(color=ng/12, linewidth=3, colorrange=(0, 1)) for ng in reverse(ngs_range)], ["$ng" for ng in reverse(ngs_range)], "Number of Groups")
    # axislegend(ax, "N Groups", nbanks=3, position=:lt)
    for filetype in ["png", "pdf"]
        save("./figures/revised/ngroups_group_level_10migrants.$filetype", fig)
    end
    fig
end
