using StaticArrays
using Random
using CairoMakie
using ProgressMeter

using Coevolution

const strategies = SVector{16,SVector{4,Bool}}(
    SA[0, 0, 0, 0], SA[1, 0, 0, 0], SA[0, 1, 0, 0], SA[1, 1, 0, 0],
    SA[0, 0, 1, 0], SA[1, 0, 1, 0], SA[0, 1, 1, 0], SA[1, 1, 1, 0],
    SA[0, 0, 0, 1], SA[1, 0, 0, 1], SA[0, 1, 0, 1], SA[1, 1, 0, 1],
    SA[0, 0, 1, 1], SA[1, 0, 1, 1], SA[0, 1, 1, 1], SA[1, 1, 1, 1],
)

begin
    Z_group = 10
    Zs = (Z_group, Z_group, Z_group, Z_group, Z_group)
    β = 1.0
    μ_s = 1 / (1 * sum(Zs))
    μ_g = 1 / (1 * sum(Zs))
    ξ = 0.9 # Likelihood of strategic update as opposed to a group update
    α = 0.9 # Assortment of interactions
    γ = 0.9 # Assortment of reproduction
    c = 1.0 # Cost of contribution
    ϵ_p = 0.01 # Error rate of production
    ϵ_c = 0.01 # Error rate of competition
    N = 2 * Z_group^2
end

S_initial = rand_S_initial_revised(Zs; strategy_set=1:16)
a = 1.0
pots = (SA[1+2., 2(1+2)], SA[1+2., 2(1+2)])
rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
@b main_simulation_loop(S_initial, N, rp; strategy_set=1:16)

l = 41
a_range = range(0, 4, length=l)
b_range = range(0, 4, length=l)

v1(b) = 1 + b
v2(b) = (1 + b)^2

# Collect only the four agnostic strategies
mean_strategy_count_matrix_grouped = let
    M = [zeros(16) for i in b_range, j in a_range]
    iterator = collect(Iterators.product(b_range, a_range))
    @showprogress Threads.@threads for ij in 1:l^2
        b, a = iterator[ij]
        S_initial = rand_S_initial_revised(Zs; strategy_set=1:16)
        pots = (SA[v1(b), v2(b)], SA[v1(b), v2(b)])
        rp = RevisedParameters(Zs, β, μ_s, μ_g, ξ, α, γ, c, a, pots, ϵ_p, ϵ_c)
        strategy_count_by_generation = main_simulation_loop(S_initial, N, rp; strategy_set=1:16)
        burn_in_period = N ÷ 10
        collection_period = N - burn_in_period
        M[ij] = dropdims(sum(strategy_count_by_generation[:, :, end-collection_period+1:end], dims=(1, 3)), dims=(1, 3)) ./ collection_period
    end
    M
end

# Debugging

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
                b_range,
                a_range,
                output_matrix[idx],
                colorrange=(0, Zs[1]),
                colormap=cmaps[idx]
            )
            # vl = vlines!(ax, [m_out], color=:black, linestyle=:dash)
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
            # save("figures/revised/α=0.7-r=[1-0.5-2mc].$filetype", fig)
        end
        display(fig)
    end
end