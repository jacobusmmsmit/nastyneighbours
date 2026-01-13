using StaticArrays
using Random
using CairoMakie

using Coevolution

strategies = SVector{16,SVector{4,Bool}}(
    SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
    SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
    SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
    SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
)

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
    cooperation_per_context_per_group_per_timestep = zeros(Int, n_groups, 8, n_timesteps)
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

begin
    Zs = (15, 15)
    β = 1
    μ = 1 / sum(Zs)
    c_p = 1
    c_c = 2.1
    m_in = 2.5
    m_out = 2.5
    α = 0.9
    ϵ_p = 0.01
    ϵ_c = 0.01
    n_migrants = maximum(Zs) ÷ 4
    inequality = 0.0
    hierarchy_strength = 0.9
    N = 4_000
end

m_range = range(1, 4, length=11)
c_range = range(0, 3, length=11)
S_initial = rand_S_initial_hierarchy(Zs)
hp = HierarchyParameters(Zs, β, μ, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c, n_migrants, inequality, hierarchy_strength)
T = main_simulation_loop(S_initial, N, hp)
@be main_simulation_loop(S_initial, N, hp)
# get_cooperation_over_time(T);

let
    fig = Figure(size=(600, 400))
    # Rows = [group 1, group 2,..., group N], columns = [out-group, in-group] 
    axs = [Axis(fig[by_group, to_relation], xlabel="Timestep", ylabel="Group $by_group", title="$(to_relation==1 ? "Out-group" : "In-group")") for (by_group, to_relation) in Iterators.product(1:length(Zs), 1:2)]
    cooperation_per_context_per_group_per_timestep = get_cooperation_over_time(T)
    for group in 1:length(Zs)
        engagement_per_context_per_timestep = @views cooperation_per_context_per_group_per_timestep[group, :, :]
        display(engagement_per_context_per_timestep)
        for (entry, is_to_ingroup, strategy) in zip(1:8, (1, 1, 1, 1, 2, 2, 2, 2), (1, 2, 3, 4, 1, 2, 3, 4))
            ax = axs[group, is_to_ingroup]
            cooperation_per_timestep = @views engagement_per_context_per_timestep[entry, :]
            lines!(ax, cooperation_per_timestep, color=Cycled(strategy))
        end
    end
    Legend(fig[1:2, 3], [LineElement(color=Cycled(i)) for i in 1:4], ["Free riders", "Producers", "Claimers", "PCs"])
    fig
end

# 4x4
S_initial = rand_S_initial_hierarchy(Zs)
mean_strategy_count_matrix_grouped = map(Iterators.product(m_range, c_range)) do (m_in, c_c)
    println("($m_in, $c_c)")
    rand_S_initial_hierarchy!(S_initial, Zs)
    hp = HierarchyParameters(Zs, β, μ, c_p, c_c, m_in, m_in, α, ϵ_p, ϵ_c, n_migrants, inequality, hierarchy_strength)
    strategy_count_by_generation = main_simulation_loop(S_initial, N, hp)
    dropdims(mean(strategy_count_by_generation, dims=(1, 3)), dims=(1, 3))
end

let
    structured_strategies = SVector{16,SVector{4,Bool}}(
        SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
        SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
        SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
        SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1]
    )


    output_matrix = reshape(
        [
            getindex.(mean_strategy_count_matrix_grouped, i)
            for i in 1:16
        ], 4, 4
    )

    cmaps = [getindex(cgrads, group) for group in [1, 2, 3, 4, 2, 2, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4]]

    begin
        figsize = (620, 600)
        fig = Figure(; size=figsize)
        gl = fig[1, 1] = GridLayout()
        axs = []
        hms = []
        for idx in 1:16
            local co, po, ci, p_i = structured_strategies[idx]

            i_claim = 5 - (1 + co + 2ci)
            j_produce = 1 + po + 2p_i
            println("$idx: ($(Int(co)), $(Int(po)), $(Int(ci)), $(Int(p_i))): at ($i_claim, $j_produce)")
            ax = Axis(gl[i_claim, j_produce]; aspect=1) #,title="$idx"
            ax.xlabel = "In-Multiplier"
            ax.ylabel = "Claiming cost"
            push!(axs, ax)
            if i_claim < 4
                ax.xticklabelsvisible = false
                ax.xticksvisible = false
                ax.xlabelvisible = false
            end
            if j_produce > 1
                ax.yticklabelsvisible = false
                ax.yticksvisible = false
                ax.ylabelvisible = false
            end
            if j_produce == 1
                label = ["Share", "Out-Claim", "In-Claim", "Uni-Claim"][i_claim]
                Label(gl[5-i_claim, 0, Makie.Right()], label; padding=(0, 5, 0, 0), rotation=π / 2, font=:bold)
            end
            if i_claim == 4
                label = ["Free-ride", "Out-Prod", "In-Prod", "Uni-Prod"][j_produce]
                Label(gl[5, j_produce, Makie.Top()], label; padding=(0, 0, 0, 5), font=:bold)
            end
            hm = heatmap!(
                ax,
                m_range,
                c_range,
                output_matrix[idx],
                colorrange=(0, Z_1 + Z_2),
                colormap=cmaps[idx]
            )
            vl = vlines!(ax, [1.5], color=:black, linestyle=:dash)
            push!(hms, hm)
        end

        Label(gl[2:3, 0, Makie.Left()], "Group-dependent Claiming", padding=(-10, 0, 0, 0), rotation=π / 2, font=:bold)
        Label(gl[5, 2:3, Makie.Bottom()], "Group-dependent Production", padding=(0, 0, -5, 0), rotation=0, font=:bold)

        Colorbar(gl[1, 5], hms[4], label="")
        Colorbar(gl[2, 5], hms[2], label="")
        Colorbar(gl[3, 5], hms[3], label="")
        Colorbar(gl[4, 5], hms[1], label="")
        Label(gl[1:4, 5, Makie.Right()], "Strategy proportion of population"; padding=(0, -50, 0, 0), rotation=3π / 2)
        yspace = maximum(tight_yticklabel_spacing!, axs)
        xspace = maximum(tight_xticklabel_spacing!, axs)
        for ax in axs
            ax.yticklabelspace = yspace
            ax.xticklabelspace = xspace
        end
        for (xbounds, ybounds) in [(1:3, 1:1), (1:3, 2:4), (4:4, 1:1), (4:4, 2:4)]
            b = Box(
                gl[xbounds, ybounds, Makie.GridLayoutBase.Outer()],
                alignmode=Outside(-7, -7, -7, -7),
                cornerradius=3,
                strokewidth=1.2,
                # linestyle=:dash,
                color=(:black, 0.0),
            )
            translate!(b.blockscene, 0, 0, -202)
        end
        highlight_horizontal = Box(
            gl[2:3, 0:4, Makie.GridLayoutBase.Outer()],
            alignmode=Outside(-10, -5, -5, -5),
            cornerradius=0,
            strokewidth=0,
            # linestyle=:dash,
            color=(:red, 0.1),
        )
        highlight_vertical = Box(
            gl[1:5, 2:3, Makie.GridLayoutBase.Outer()],
            alignmode=Outside(-3, -3, -10, 20),
            cornerradius=0,
            strokewidth=0,
            # linestyle=:dash,
            color=(:red, 0.1),
        )
        translate!(highlight_horizontal.blockscene, 0, 0, -200)
        translate!(highlight_vertical.blockscene, 0, 0, -201)
        # highlight_box = Box(
        #     gl[2:3, 1:4, Makie.GridLayoutBase.Outer()], alignmode=Outside(-10, -10, -10, -10, cornerradius=3, strokewidth=1, color=(:red, 0.1))
        # )
        for i in 1:4
            # colsize!(gl, i, Relative(0.16))
            # rowsize!(gl, i, Relative(0.19))
        end
        # colsize!(gl, 5, Aspect(1, 0.3))
        local small = Relative(0.03)
        local big = Relative(0.07)
        colgap!(gl, 1, small)
        colgap!(gl, 2, big)
        colgap!(gl, 3, small)
        colgap!(gl, 4, small)
        rowgap!(gl, 1, small)
        rowgap!(gl, 2, small)
        rowgap!(gl, 3, big)
        colsize!(gl, 0, Relative(0.05))
        rowsize!(gl, 5, Relative(0.05))
        label_options = (;
            padding=(0, 0, 10, 0),
            justification=:left,
            halign=:left,
            font=:bold,
        )
        Label(gl[1, 1, Makie.Top()], "A: Claimers"; label_options..., padding=(-35, 0, 10, 0))
        Label(gl[1, 2:4, Makie.Top()], "B: Produce-Claimers"; label_options...)
        Label(gl[4, 1, Makie.Top()], "C: Freeriders"; label_options..., padding=(-35, 0, 10, 0))
        Label(gl[4, 2:4, Makie.Top()], "D: Producers"; label_options...)
        # Colorbar(fig[:, 3], hms[1], colorrange=(0, 1), label="Number of agents")
        for filetype in ("png", "pdf")
            save("figures/hierarchy/inequality_$(inequality)_$hierarchy_strength.$filetype", fig)
        end
        display(fig)
    end
end