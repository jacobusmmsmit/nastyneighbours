using StaticArrays
using Random
using CairoMakie

includet("../implementation/common.jl")
includet("../implementation/payoff_from_interaction.jl")
includet("../implementation/structured.jl")

begin
    Z_1 = 100
    Z_2 = 0
    β = 1
    μ = 1 / Z_1 + Z_2
    c_p = 1
    c_c = 1
    m_in = 10000
    m_out = 1.5
    α = 0.5
    ϵ_p = 0
    ϵ_c = 0
    N = 100#_000
end

m_range = range(1, 4, length=31)
c_range = range(0, 3, length=31)
S_initial = rand_S_initial_structured(Z_1, Z_2)
sp = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c)

# strategy_count_by_generation = main_simulation_loop(S_initial, N, sp)

# let
#     fig = Figure(size=(600, 400))
#     ax = Axis(fig[1, 1])
#     ax2 = Axis(fig[2, 1])
#     local m_e = 1.5
#     local c_c_e = 0.2
#     S_initial = rand_S_initial_structured(Z_1, Z_2)
#     sp_e = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c_e, m_e, m_e, 1, ϵ_p, ϵ_c)
#     strategy_count_by_generation_e = main_simulation_loop(S_initial, 1000, sp_e)
#     for row_i in 0:3
#         row = vec(sum(strategy_count_by_generation_e[(1:4).+4row_i, :], dims=1))
#         lines!(ax, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
#     end
#     local m_f = 2.5
#     local c_c_f = 2.1
#     sp_f = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c_f, m_f, m_f, 1, ϵ_p, ϵ_c)
#     strategy_count_by_generation_f = main_simulation_loop(S_initial, 1000, sp_f)
#     for row_i in 0:3
#         row = vec(sum(strategy_count_by_generation_f[(1:4).+4row_i, :], dims=1))
#         lines!(ax2, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
#     end
#     display(fig)
# end

# mean_strategy_count_matrix = map(Iterators.product(m_range, c_range)) do (m, c_c)
#     println("($m, $c_c)")
#     S_initial = rand_S_initial_structured(Z_1, Z_2)
#     sp = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c, m, m, α, ϵ_p, ϵ_c)
#     strategy_count_by_generation = main_simulation_loop(S_initial, N, sp)
#     [sum(mean.(rows)) for rows in Iterators.partition(eachrow(strategy_count_by_generation), 4)][1:4]
# end

# begin
#     figsize = (1100, 600)
#     fig2 = Figure(; size=figsize)
#     ga = fig2[1, 1] = GridLayout()
#     gb = fig2[1, 2] = GridLayout()
#     begin
#         axs2a = []
#         hms2a = []
#         for i in 1:2, j in 1:2
#             idx = 2(i - 1) + (j - 1)
#             ax = Axis(ga[3-j, 2i-1]; aspect=1)
#             push!(axs2a, ax)
#             ax.xlabel = "Productivity"
#             ax.ylabel = "Claiming cost"
#             ax.title = labels[idx+1]
#             hm = heatmap!(
#                 ax,
#                 m_range,
#                 c_range,
#                 getindex.(mean_strategy_count_matrix, idx + 1),
#                 colorrange=(0, Z_1 + Z_2),
#                 colormap=cgrads[idx+1],
#             )
#             push!(hms2a, hm)
#             cb = Colorbar(ga[3-j, 2i], hm, label="Number of agents", tellheight=true)
#             cb.height = Relative(0.73)
#         end
#         for (label, pos) in zip(["a", "b", "c", "d"], [[1, 1], [1, 3], [2, 1], [2, 3]])
#             Label(ga[pos[1], pos[2], TopLeft()], label,
#                 fontsize=26,
#                 font=:bold,
#                 padding=(0, 5, 0, 0),
#                 halign=:right)
#         end
#     end
#     ax = Axis(gb[1, 1])
#     ax2 = Axis(gb[2, 1])
#     local m_e = 1.5
#     local c_c_e = 0.2
#     S_initial = rand_S_initial_structured(Z_1, Z_2)
#     sp_e = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c_e, m_e, m_e, 1, ϵ_p, ϵ_c)
#     strategy_count_by_generation_e = main_simulation_loop(S_initial, 1000, sp_e)
#     for row_i in 0:3
#         row = vec(sum(strategy_count_by_generation_e[(1:4).+4row_i, :], dims=1))
#         lines!(ax, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
#     end
#     local m_f = 2.5
#     local c_c_f = 2.1
#     sp_f = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c_f, m_f, m_f, 1, ϵ_p, ϵ_c)
#     strategy_count_by_generation_f = main_simulation_loop(S_initial, 1000, sp_f)
#     for row_i in 0:3
#         row = vec(sum(strategy_count_by_generation_f[(1:4).+4row_i, :], dims=1))
#         lines!(ax2, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
#     end
#     # display(fig)
#     elements = [MarkerElement(; marker=:rect, color=color, markersize=20) for color in strat_colours]
#     Legend(gb[3, 1], elements, labels, "Strategies", orientation=:horizontal)
#     Label(gb[1, 1, TopLeft()], "e",
#             fontsize=26,
#             font=:bold,
#             padding=(0, 15, 5, 0),
#             halign=:right)
#     Label(gb[2, 1, TopLeft()], "f",
#         fontsize=26,
#         font=:bold,
#         padding=(0, 15, 10, 0),
#         halign=:right)
#     fax = Axis(ga[3, :])
#     hidedecorations!(fax)  # hides ticks, grid and lables
#     hidespines!(fax)
#     rowsize!(ga, 1, Relative(0.42))
#     rowsize!(ga, 2, Relative(0.42))
#     rowgap!(ga, 0)
#     for filetype in ("png", "pdf")
#         save("figures/fig2.$filetype", fig2)
#     end
#     display(fig2)
# end

for α in 0.5:0.1:0.9
    # 4x4
    Z_1 = 50
    Z_2 = 50
    m_out = 1.5
    N = 100_000
    S_initial_grouped = rand_S_initial_structured(Z_1, Z_2)
    mean_strategy_count_matrix_grouped = map(Iterators.product(m_range, c_range)) do (m_in, c_c)
        println("($m_in, $c_c)")
        rand_S_initial_structured!(S_initial, Z_1, Z_2)
        sp = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c)
        strategy_count_by_generation = main_simulation_loop(S_initial, N, sp)
        v = vec(mean(strategy_count_by_generation; dims=2))
        sum(v[1:16]) ≈ sum(v[17:32]) || error("$(sum(v[1:16])) != $(sum(v[17:32]))")
        v[1:16] + v[17:32]
    end




    # sj_claim_out, sj_produce_out, sj_claim_in, sj_produce_in
    structured_strategies = SVector{16,SVector{4,Bool}}(
        SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
        SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
        SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
        SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1]
    )


    output_matrix = reshape([
            getindex.(mean_strategy_count_matrix_grouped, i)
            for i in 1:16
        ], 4, 4)

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
            ax = Axis(gl[i_claim, j_produce]; aspect=1)#,title="$idx"
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
            save("figures/N4b4_$α.$filetype", fig)
        end
        display(fig)
    end
end