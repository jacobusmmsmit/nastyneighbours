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
    N = 2 * Z_group^2
end

S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
a = 1.0
pots = (SA[1+2., 2(1+2)], SA[1+2., 2(1+2)])
rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
@b main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)

l = 41
a_range = range(0, 4, length=l)
b_range = range(0, 4, length=l)

v1(b) = 1 + b
v2(b) = 2(1 + b)

# Collect only the four agnostic strategies
mean_strategy_count_matrix_grouped = let
    M = [zeros(4) for i in b_range, j in a_range]
    iterator = collect(Iterators.product(b_range, a_range))
    @showprogress Threads.@threads for ij in 1:l^2
        b, a = iterator[ij]
        S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
        pots = (SA[v1(b), v2(b)], SA[v1(b), v2(b)])
        rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
        strategy_count_by_generation = main_simulation_loop(S_initial, N, rp; strategy_set=group_agnostic_strategies)
        burn_in_period = N ÷ 10
        collection_period = N - burn_in_period
        M[ij] = dropdims(sum(strategy_count_by_generation[:, :, end-collection_period+1:end], dims=(1, 3)), dims=(1, 3))[[1, 6, 11, 16]] ./ collection_period
    end
    M
end

begin
    fig2 = Figure(; size=(1100, 600))
    ga = fig2[1, 1] = GridLayout()
    gb = fig2[1, 2] = GridLayout()
    begin
        for i in 1:2, j in 1:2
            idx = 2(i - 1) + (j - 1)
            ax = Axis(ga[3-j, 2i-1]; aspect=1)
            ax.xlabel = "Productivity"
            ax.ylabel = "Claiming cost"
            ax.title = labels[idx+1]
            hm = heatmap!(
                ax,
                b_range,
                a_range,
                getindex.(mean_strategy_count_matrix_grouped, idx + 1),
                colorrange=(0, sum(Zs)),
                colormap=cgrads[idx+1],
            )
            cb = Colorbar(ga[3-j, 2i], hm, label="Number of agents", tellheight=true)
            cb.height = Relative(0.73)
            # Freeriders polygon
            poly!(
                [Point2f(-1, 1 / 2), Point2f(3 / 2, 1 / 2), Point2f(2, 1), Point2f(2, 10), Point2f(-1, 10)],
                strokecolor=(strat_colours[1], 1),
                strokewidth=2,
                # linestyle=:dash,
                color=(:black,0.0)
            )
            # Claiming polygon
            poly!(
                [Point2f(-1, -1), Point2f(2, -1), Point2f(2, 1 / 2), Point2f(-1, 1 / 2)],
                strokecolor=(strat_colours[2], 1),
                strokewidth=2,
                # linestyle=:dash,
                color=(:black, 0.0),
            )
            # PC
            poly!(
                [Point2f(1, -1), Point2f(1, 1), Point2f(6, 6), Point2f(6, -1)],
                strokecolor=(strat_colours[4], 1),
                strokewidth=2,
                # linestyle=:dash,
                color=(:black, 0.0),
            )
            # Regions for v1 = 1+b, v2 = (1+b)^2
            # poly!(
            #     [Point2f(-1, 1 / 2), Point2f(3 / 2, 1 / 2), Point2f(2, 1), Point2f(2, 10), Point2f(-1, 10)],
            #     strokecolor=(strat_colours[1], 1),
            #     strokewidth=2,
            #     # linestyle=:dash,
            #     color=(:black,  0.0)
            # )
            # poly!(
            #     [Point2f(-1, -1), Point2f(2, -1), Point2f(2, 1 / 2), Point2f(-1, 1 / 2)],
            #     strokecolor=(strat_colours[2], 1),
            #     strokewidth=2,
            #     # linestyle=:dash,
            #     color=(:black, 0.0),
            # )
            # poly!(
            #     [Point2f(1, -3), Point2f.(1:0.25:5, (b -> (1 / 2) * (b^2 + 2b - 1)).(1:0.25:5))..., Point2f(5, -3)],
            #     strokecolor=(strat_colours[4], 1),
            #     strokewidth=2,
            #     # linestyle=:dash,
            #     color=(:black, 0.0),
            # )
            limits!(ax, (0, 4), (0, 4))
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
    # error()
    ax = Axis(gb[1, 1])
    ax2 = Axis(gb[2, 1])
    b_e = 1.2
    a_e = 0.7
    N_e = Z_group^2
    S_initial = rand_S_initial_revised(Zs; strategy_set=group_agnostic_strategies)
    pots_e = (SA[v1(b_e), v2(b_e)], SA[v1(b_e), v2(b_e)])
    rp_e = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots_e, ϵ_p, ϵ_c)
    strategy_count_by_generation_e = main_simulation_loop(S_initial, N_e, rp_e; strategy_set=group_agnostic_strategies)
    for row_i in 0:3
        row = vec(sum(strategy_count_by_generation_e[:, (1:4).+4row_i, :], dims=(1, 2)))
        lines!(ax, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
    end
    b_f = 2.5
    a_f = 3.5
    pots_f = (SA[v1(b_f), v2(b_f)], SA[v1(b_f), v2(b_f)])
    rp_f = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a_f, pots_f, ϵ_p, ϵ_c)
    strategy_count_by_generation_f = main_simulation_loop(S_initial, N_e, rp_f; strategy_set=group_agnostic_strategies)
    for row_i in 0:3
        row = vec(sum(strategy_count_by_generation_f[:, (1:4).+4row_i, :], dims=(1, 2)))
        lines!(ax2, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
    end
    elements = [MarkerElement(; marker=:rect, color=color, markersize=20) for color in strat_colours]
    Legend(gb[3, 1], elements, labels, "Strategies", orientation=:horizontal)
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
    fax = Axis(ga[3, :])
    hidedecorations!(fax)  # hides ticks, grid and lables
    hidespines!(fax)
    rowsize!(ga, 1, Relative(0.42))
    rowsize!(ga, 2, Relative(0.42))
    rowgap!(ga, 0)
    for filetype in ("png", "pdf")
        save("figures/revised/fig2_v1_1+b_v2_2(1+b)_$(N)_$(l).$filetype", fig2)
    end
    display(fig2)
end



# cmaps = [getindex(cgrads, group) for group in [1, 2, 3, 4, 2, 2, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4]]

# begin
#     figsize = (620, 600)
#     fig = Figure(; size=figsize)
#     gl = fig[1, 1] = GridLayout()
#     axs = []
#     hms = []
#     for idx in 1:16
#         local co, po, ci, p_i = strategies[idx]
#         i_claim = 5 - (1 + co + 2ci)
#         j_produce = 1 + po + 2p_i
#         # println("$idx: ($(Int(co)), $(Int(po)), $(Int(ci)), $(Int(p_i))): at ($i_claim, $j_produce)")
#         ax = Axis(gl[i_claim, j_produce]; aspect=1) #,title="$idx"
#         ax.xlabel = "In-Multiplier"
#         ax.ylabel = "Claiming cost"
#         push!(axs, ax)
#         if i_claim < 4
#             ax.xticklabelsvisible = false
#             ax.xticksvisible = false
#             ax.xlabelvisible = false
#         end
#         if j_produce > 1
#             ax.yticklabelsvisible = false
#             ax.yticksvisible = false
#             ax.ylabelvisible = false
#         end
#         if j_produce == 1
#             label = ["Share", "Out-Claim", "In-Claim", "Uni-Claim"][i_claim]
#             Label(gl[5-i_claim, 0, Makie.Right()], label; padding=(0, 5, 0, 0), rotation=π / 2, font=:bold)
#         end
#         if i_claim == 4
#             label = ["Free-ride", "Out-Prod", "In-Prod", "Uni-Prod"][j_produce]
#             Label(gl[5, j_produce, Makie.Top()], label; padding=(0, 0, 0, 5), font=:bold)
#         end
#         hm = heatmap!(
#             ax,
#             b_range,
#             a_range,
#             getindex.(mean_strategy_count_matrix_grouped, idx),
#             colorrange = (0, sum(Zs)),
#             colormap=cmaps[idx]
#         )
#         vl = vlines!(ax, [1.5], color=:black, linestyle=:dash)
#         push!(hms, hm)
#     end

#     Label(gl[2:3, 0, Makie.Left()], "Group-dependent Claiming", padding=(-10, 0, 0, 0), rotation=π / 2, font=:bold)
#     Label(gl[5, 2:3, Makie.Bottom()], "Group-dependent Production", padding=(0, 0, -5, 0), rotation=0, font=:bold)

#     Colorbar(gl[1, 5], hms[4], label="")
#     Colorbar(gl[2, 5], hms[2], label="")
#     Colorbar(gl[3, 5], hms[3], label="")
#     Colorbar(gl[4, 5], hms[1], label="")
#     Label(gl[1:4, 5, Makie.Right()], "Strategy proportion of population"; padding=(0, -50, 0, 0), rotation=3π / 2)
#     yspace = maximum(tight_yticklabel_spacing!, axs)
#     xspace = maximum(tight_xticklabel_spacing!, axs)
#     for ax in axs
#         ax.yticklabelspace = yspace
#         ax.xticklabelspace = xspace
#     end
#     for (xbounds, ybounds) in [(1:3, 1:1), (1:3, 2:4), (4:4, 1:1), (4:4, 2:4)]
#         b = Box(
#             gl[xbounds, ybounds, Makie.GridLayoutBase.Outer()],
#             alignmode=Outside(-7, -7, -7, -7),
#             cornerradius=3,
#             strokewidth=1.2,
#             # linestyle=:dash,
#             color=(:black, 0.0),
#         )
#         translate!(b.blockscene, 0, 0, -202)
#     end
#     highlight_horizontal = Box(
#         gl[2:3, 0:4, Makie.GridLayoutBase.Outer()],
#         alignmode=Outside(-10, -5, -5, -5),
#         cornerradius=0,
#         strokewidth=0,
#         # linestyle=:dash,
#         color=(:red, 0.1),
#     )
#     highlight_vertical = Box(
#         gl[1:5, 2:3, Makie.GridLayoutBase.Outer()],
#         alignmode=Outside(-3, -3, -10, 20),
#         cornerradius=0,
#         strokewidth=0,
#         # linestyle=:dash,
#         color=(:red, 0.1),
#     )
#     translate!(highlight_horizontal.blockscene, 0, 0, -200)
#     translate!(highlight_vertical.blockscene, 0, 0, -201)
#     # highlight_box = Box(
#     #     gl[2:3, 1:4, Makie.GridLayoutBase.Outer()], alignmode=Outside(-10, -10, -10, -10, cornerradius=3, strokewidth=1, color=(:red, 0.1))
#     # )
#     for i in 1:4
#         # colsize!(gl, i, Relative(0.16))
#         # rowsize!(gl, i, Relative(0.19))
#     end
#     # colsize!(gl, 5, Aspect(1, 0.3))
#     local small = Relative(0.03)
#     local big = Relative(0.07)
#     colgap!(gl, 1, small)
#     colgap!(gl, 2, big)
#     colgap!(gl, 3, small)
#     colgap!(gl, 4, small)
#     rowgap!(gl, 1, small)
#     rowgap!(gl, 2, small)
#     rowgap!(gl, 3, big)
#     colsize!(gl, 0, Relative(0.05))
#     rowsize!(gl, 5, Relative(0.05))
#     label_options = (;
#         padding=(0, 0, 10, 0),
#         justification=:left,
#         halign=:left,
#         font=:bold,
#     )
#     Label(gl[1, 1, Makie.Top()], "A: Claimers"; label_options..., padding=(-35, 0, 10, 0))
#     Label(gl[1, 2:4, Makie.Top()], "B: Produce-Claimers"; label_options...)
#     Label(gl[4, 1, Makie.Top()], "C: Freeriders"; label_options..., padding=(-35, 0, 10, 0))
#     Label(gl[4, 2:4, Makie.Top()], "D: Producers"; label_options...)
#     # Colorbar(fig[:, 3], hms[1], colorrange=(0, 1), label="Number of agents")
#     for filetype in ("png", "pdf")
#         # save("figures/N4b4_$α.$filetype", fig)
#     end
#     display(fig)
# end


let
    begin
        Z_group = 30
        Zs = (Z_group,)
        β = 1.0
        μ_s = 1 / (10 * sum(Zs))
        μ_g = 1 / (1 * sum(Zs))
        ξ = 1.0 # Likelihood of strategic update as opposed to a group update
        α = 1.0 # Assortment of interactions
        γ = 1.0 # Assortment of reproduction
        c = 1.0 # Cost of contribution
        ϵ_p = 0.01 # Error rate of production
        ϵ_c = 0.01 # Error rate of competition
        # N = 50
    end
    S_initial = rand_S_initial_revised(Zs; strategy_set=SA[1, 6, 16])
    b = 1.3
    a = 0.5
    pots = (SA[1+b, 2(1+b)], SA[1+b, 2(1+b)])
    rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
    strategy_count_by_generation = main_simulation_loop(S_initial, N, rp; strategy_set=SA[1, 6, 16])
    fig = Figure()
    ax = Axis(fig[1, 1])
    for row_i in 0:3
        row = vec(sum(strategy_count_by_generation[:, (1:4).+4row_i, :], dims=(1, 2)))
        lines!(ax, row, linewidth=3, alpha=1, label="$(row_i - 16)", color=strat_colours[row_i+1])
    end
    fig
end

let
    begin
        Z_group = 30
        Zs = (Z_group,)
        β = 1.0
        μ_s = 1 / (10 * sum(Zs))
        μ_g = 1 / (1 * sum(Zs))
        ξ = 1.0 # Likelihood of strategic update as opposed to a group update
        α = 1.0 # Assortment of interactions
        γ = 1.0 # Assortment of reproduction
        c = 1.0 # Cost of contribution
        ϵ_p = 0.01 # Error rate of production
        ϵ_c = 0.01 # Error rate of competition
        # N = 50
    end
    iterator = [(i, j) for i in 0:Z_group for j in 0:Z_group-i]
    for (n_freeriders, n_produceclaimers) in iterator
        n_claimers = Z_group - (n_freeriders + n_produceclaimers)
        S_initial = reshape([n_freeriders, 0, 0, 0, 0, n_claimers, 0, 0, 0, 0, 0, 0, 0, 0, 0, n_produceclaimers], 1, 16)
        n_freeriders_end, n_produceclaimers_end = main_simulation_loop(S_initial, 10, rp; strategy_set=SA[1, 6, 16])[1, [1, 16], end]
    end
end

let
    begin
        # Z_group = 
        Zs = (Z_group,)
        β = 1.0
        μ_s = 1 / (1000 * sum(Zs))
        μ_g = 1 / (1 * sum(Zs))
        ξ = 1.0 # Likelihood of strategic update as opposed to a group update
        α = 1.0 # Assortment of interactions
        γ = 1.0 # Assortment of reproduction
        c = 1.0 # Cost of contribution
        ϵ_p = 0.01 # Error rate of production
        ϵ_c = 0.01 # Error rate of competition
        N = 20
    end
    b_e = 1.2
    a_e = 0.7
    b_f = 2.5
    a_f = 3.5
    local b = b_f
    local a = a_f
    iterator = [(i, j) for i in 0:Z_group for j in 0:Z_group-i]
    dirs = let
        M = zeros(Point2d, Z_group + 1, Z_group + 1)
        @showprogress Threads.@threads for (n_freeriders, n_produceclaimers) in iterator
            from = SA[n_freeriders, n_produceclaimers]
            n_claimers = Z_group - (n_freeriders + n_produceclaimers)
            S_initial = reshape([n_freeriders, 0, 0, 0, 0, 0, 0, 0, 0, 0, n_claimers, 0, 0, 0, 0, n_produceclaimers], 1, 16)
            pots = (SA[1+b, 2(1+b)], SA[1+b, 2(1+b)])
            rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
            n_freeriders_end, n_produceclaimers_end = mean(main_simulation_loop(S_initial, N, rp; strategy_set=SA[1, 11, 16])[1, [1, 16], end] for _ in 1:N^2)
            dir = Point2d(n_freeriders_end - n_freeriders, n_produceclaimers_end - n_produceclaimers)
            M[n_freeriders+1, n_produceclaimers+1] = dir
        end
        M
    end

    interp_linear = linear_interpolation(
        (0:Z_group, 0:Z_group),
        dirs / maximum(norm.(dirs)),
        extrapolation_bc=SA[0, 0]
    )
    # function f1(x, y, Z_group)
    #     shear = SA[1 0.5; 0 1]
    #     # Transform back to (n_freeriders, n_produceclaimers)-space
    #     original_x, original_y = inv(shear) * SA[x, y]
    #     # Perform the interpolation
    #     original_dir = interp_linear(original_x, original_y)
    #     # Transform it to the ternary space
    #     return Point2f(shear * original_dir)
    # end
    fig = Figure(size=(500, 500))
    ax = Axis(
        fig[1, 1],
        aspect=1,
        title="Dynamics of Freeriders, Producers, and Produce-Claimers",
        subtitle="Benefit: $b, Aggression cost: $a",
        xlabel="Freeriders",
        ylabel="ProduceClaimers"
    )
    # hidedecorations!(ax)
    # hidespines!(ax)
    sp = streamplot!(
        ax,
        ((x, y)) -> interp_linear(x, y),
        0 .. Z_group,
        0 .. Z_group,
        # colormap=streamplot_colormap,
        # colorrange=(0, 1),
        density=0.5,
    )
    scatter!([Point2f(0, Z_group), Point2f(Z_group, 0), Point2f(0, 0)], color=strat_colours[[4, 1, 3]], marker=:circle, markersize=25)
    # linesegments!(
    #     ax,
    #     Z_group .* (1 / 2Z_group .+ [0, 0.5, 0.5, 1, 1, 0]),
    #     Z_group .* [0, 1, 1, 0, 0, 0],
    #     color=:black
    # )
    # text!(
    #     ax,
    #     [Z_group .* Point2f(0.5, -0.075), Z_group .* Point2f(0.25 - 0.0375, 0.5 + 0.0375), Z_group .* Point2f(0.75 + 0.075, 0.5 + 0.075)],
    #     text=["Number of Freeriders", "Number of ProduceClaimers", "Number of Claimers"],
    #     rotation=[0, π / 3, -π / 3],
    #     align=(:center, :center)
    # )
    save("figures/revised/streamplot_b$(b)_a$(a).png", fig)
    fig
end



# let
#     hear = SA[1 0.5; 0 1]
#     (x, y) = [25, 10]
#     # Transform back to (n_freeriders, n_produceclaimers)-space
#     original_x, original_y = inv(shear) * SA[Z_group - x, y]
#     # Perform the interpolation
#     # Transform it to the ternary space
#     untranslated_x, untranslated_y = shear * SA[original_x, original_y]
#     Point2f(Z_group - untranslated_x, untranslated_y)
# end


# Testing
let
    b = 2.5
    a = 2.5
    pot = SA[1+b, 2(1+b)]
    rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
    get_payoff_matrix(pot, c, a, 1 // 2)
end