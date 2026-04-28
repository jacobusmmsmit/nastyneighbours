using CairoMakie
using Coevolution: strat_colours

begin
    fig = Figure(size=(400, 430))
    ax = Axis(
        fig[1, 1],
        aspect=1,
        xlabel="Benefit (b)",
        ylabel="Cost of Aggression (a)",
        # title="Regions where each strategy is ESS",
        subtitle="Parametrisation: v₁ = 1+b, v₂ = (1+b)²",
        titlealign=:left,
    )

    xlims!(ax, 0, 5)
    ylims!(ax, 0, 5)
    poly!(
        [Point2f(-1, 1 / 2), Point2f(3 / 2, 1 / 2), Point2f(2, 1), Point2f(2, 10), Point2f(-1, 10)],
        strokecolor=:black,
        strokewidth=2,
        linestyle=:dash,
        alpha=0.5,
        color=strat_colours[1]
    )
    poly!(
        [Point2f(-1, -1), Point2f(2, -1), Point2f(2, 1 / 2), Point2f(-1, 1 / 2)],
        strokecolor=:black,
        strokewidth=2,
        linestyle=:dash,
        alpha=0.5,
        color=strat_colours[2],
    )
    poly!(
        [Point2f(1, -3), Point2f.(1:0.25:5, (b -> (1 / 2) * (b^2 + 2b - 1)).(1:0.25:5))..., Point2f(5, -3)],
        strokecolor=:black,
        strokewidth=2,
        linestyle=:dash,
        alpha=0.5,
        color=strat_colours[4],
    )
    text!(ax, [(1.0, 3.0), (0.85, 0.25), (3.5, 2.5)], text=["N", "C", "PC"], align=(:center, :center))
    for filetype in ["pdf", "png"]
        save("figures/revised/ESS_regions_v2=(1+b)^2.$filetype", fig)
    end
    fig
end

begin
    fig = Figure(size=(400, 430))
    ax = Axis(
        fig[1, 1],
        aspect=1,
        xlabel="Benefit (b)",
        ylabel="Cost of Aggression (a)",
        # title="Regions where each strategy is ESS",
        subtitle="Parametrisation: v₁ = 1+b, v₂ = 2(1+b)",
        titlealign=:left,
    )

    xlims!(ax, 0, 5)
    ylims!(ax, 0, 5)
    # N
    poly!(
        [Point2f(-1, 1 / 2), Point2f(3 / 2, 1 / 2), Point2f(2, 1), Point2f(2, 10), Point2f(-1, 10)],
        strokecolor=:black,
        strokewidth=2,
        linestyle=:dash,
        alpha=0.5,
        color=strat_colours[1]
    )
    # C
    poly!(
        [Point2f(-1, -1), Point2f(2, -1), Point2f(2, 1 / 2), Point2f(-1, 1 / 2)],
        strokecolor=:black,
        strokewidth=2,
        linestyle=:dash,
        alpha=0.5,
        color=strat_colours[2],
    )
    # PC
    poly!(
        [Point2f(1, -1), Point2f(1, 1), Point2f(6, 6), Point2f(6, -1)],
        strokecolor=:black,
        strokewidth=2,
        linestyle=:dash,
        alpha=0.5,
        color=strat_colours[4],
    )
    text!(ax, [(1.0, 3.0), (0.85, 0.25), (3.5, 2.5)], text=["N", "C", "PC"], align=(:center, :center))
    for filetype in ["pdf", "png"]
        save("figures/revised/ESS_regions_v2=2+2b.$filetype", fig)
    end
    fig
end