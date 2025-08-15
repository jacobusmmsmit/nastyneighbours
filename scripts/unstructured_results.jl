using StaticArrays
using Random
using CairoMakie

include("../implementation/common.jl")
include("../implementation/payoff_from_interaction.jl")
include("../implementation/unstructured.jl")



begin
    Z = 100
    β = 1
    μ = 1 / Z
    c_p = 1
    c_c = 1
    m = 2
    ϵ_p = 0
    ϵ_c = 0
    S_initial = rand_S_initial_unstructured(Z)
    up = UnstructuredParameters(Z, β, μ, c_p, c_c, m, ϵ_p, ϵ_c)
    N = 1_0
end

m_range = range(1, 4, length=31)
c_range = range(0, 3, length=31)
S_initial = rand_S_initial_unstructured(Z)

mean_strategy_count_matrix = map(Iterators.product(m_range, c_range)) do (m, c_c)
    println("($m, $c_c)")
    rand_S_initial_unstructured!(S_initial, Z)
    up = UnstructuredParameters(Z, β, μ, c_p, c_c, m, ϵ_p, ϵ_c)
    strategy_count_by_generation = main_simulation_loop(S_initial, N, up)
    vec(mean(strategy_count_by_generation; dims=2))
end

begin
    begin
        figsize = (1100, 600)
        fig2 = Figure(; size=figsize)
        ga = fig2[1, 1] = GridLayout()
        gb = fig2[1, 2] = GridLayout()
        begin
            axs2a = []
            hms2a = []
            for i in 1:2, j in 1:2
                idx = (i - 1) + 2(j - 1)
                ax = Axis(ga[3-i, 2j-1]; aspect=1)
                push!(axs2a, ax)
                ax.xlabel = "Productivity"
                ax.ylabel = "Claiming cost"
                ax.title = labels[idx+1]
                hm = heatmap!(
                    ax,
                    m_range,
                    c_range,
                    getindex.(mean_strategy_count_matrix, idx + 1),
                    colorrange=(0, Z),
                    colormap=cgrads[idx+1],
                )
                push!(hms2a, hm)
                cb = Colorbar(ga[3-i, 2j], hm, label="Number of agents", tellheight=true)
                cb.height = Relative(0.73)
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
    N_lineplots = 1_000
    begin
        ax = Axis(gb[2, 1], title="Productivity: 2.5, Claiming Cost: 2.1", xlabel="Generation", ylabel="Strategy prevalence")
        m_e = 2.5
        c_c_e = 2.1
        up_e = UnstructuredParameters(Z, β, μ, c_p, c_c_e, m_e, ϵ_p, ϵ_c)
        S_initial = rand_S_initial_unstructured(Z)
        for (row_i, row) in enumerate(eachrow(main_simulation_loop(S_initial, N_lineplots, up_e)))
            lines!(ax, 1:N_lineplots, row, color=strat_colours[row_i], linewidth=3, alpha=1)
        end
    end
    begin
        ax = Axis(gb[1, 1], title="Productivity: 1.5, Claiming Cost: 0.2", xlabel="Generation", ylabel="Strategy prevalence")
        m_f = 1.5
        c_c_f = 0.2
        up_f = UnstructuredParameters(Z, β, μ, c_p, c_c_f, m_f, ϵ_p, ϵ_c)
        S_initial = rand_S_initial_unstructured(Z)
        for (row_i, row) in enumerate(eachrow(main_simulation_loop(S_initial, N_lineplots, up_f)))
            lines!(ax, 1:N_lineplots, row, color=strat_colours[row_i], linewidth=3, alpha=1)
        end
    end
    Label(gb[1, 1, TopLeft()], "e",
        fontsize=26,
        font=:bold,
        padding=(0, 15, 5, 0),
        halign=:right)
    Label(gb[2, 1, TopLeft()], "f",
        fontsize=26,
        font=:bold,
        padding=(0, 15, 10, 0),
        halign=:right)
    fax = Axis(ga[3, 1:2])
    hidedecorations!(fax)  # hides ticks, grid and lables
    hidespines!(fax)
    rowsize!(ga, 1, Relative(0.42))
    rowsize!(ga, 2, Relative(0.42))
    rowgap!(ga, 0)
    elements = [MarkerElement(; marker=:rect, color=color, markersize=20) for color in strat_colours]
    Legend(gb[3, 1], elements, labels, "Strategies", orientation=:horizontal)
    for filetype in ("png", "pdf")
        save("figures/fig2.$filetype", fig2)
    end
    display(fig2)
end