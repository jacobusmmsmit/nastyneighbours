using CairoMakie

begin
    figsize = (550, 450)
    fig = Figure(; size=figsize)
    ax = Axis(fig[1, 1];
        aspect=1,
        title="Regions where each strategy is an ESS",
        limits = ((1,4), (0, 3))
        # subtitle="(Simulated vs Analytical)"
    )
    ax.xlabel = "Contribution multiplier"
    ax.ylabel = "Claiming cost"
    poly!(ax, Point2f[(1, 0), (2, 1), (2, 3), (1, 3)], color=(:black, 0.5),strokecolor=:black,linestyle=:dash,strokewidth=2,label=labels[1])
    poly!(ax, Point2f[(2, 2), (3, 3), (2, 3)], color=(strat_colours[3], 0.5),strokecolor=strat_colours[3],linestyle=:dash,strokewidth=2,label=labels[3])
    poly!(ax, Point2f[(2, 0), (2, 1), (4, 3), (4, 0)], color=(strat_colours[4], 0.5),strokecolor=strat_colours[4],linestyle=:dash,strokewidth=2,label=labels[4])
    elements = [
        PolyElement(
            color= (fill_colour, 0.5),
            strokecolor= fill_colour, strokewidth=2,linestyle=:dash
            )
        for fill_colour in strat_colours
    ]
    # annotation!(ax, [(1.5, 2), ], text= [L"m < \min\{2, 1+\frac{c_c}{c_p}\}", ], color=:white)
    Legend(fig[1,2], elements, labels)
    for filetype in ("png", "pdf")
        save("figures/ESS_regions.$filetype", fig)
    end
    display(fig)

end