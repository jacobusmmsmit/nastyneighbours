using StaticArrays
using CairoMakie

include("../implementation/common.jl")
include("../implementation/structured.jl")

one_at(i, n) = (x = zeros(n); x[i] = 1.0; return x)
to_bin(i) = SVector{2,Bool}(i >> shift & 1 != 0 for shift in false:true)

function nash_equilibrium(m, c_c, c_p)
    c = c_c / c_p
    if (m < 2)
        if (m < c + 1)
            return one_at(1, 4)
        else

            return SA[2/m-1, 2-2(c+1)/m, 0, 2c/m]
        end
    else
        if m < c
            return one_at(3, 4)
        elseif m > c + 1
            return one_at(4, 4)
        else
            return SA[2(1-c/m), 0, 2(1+c-m)/m, 1-2/m]
        end
    end
end


begin
    figsize = (550, 450)
    fig = Figure(; size=figsize)
    axs = map(1:4) do idx
        _i, _j = to_bin(idx - 1)
        i = 2 - _i
        j = 2 * _j + 1
        ax = Axis(fig[i, j];
            aspect=1,
            title=labels[idx],
            # subtitle="(Simulated vs Analytical)"
        )
        ax.xlabel = "Contribution multiplier"
        ax.ylabel = "Claiming cost"
        ax
    end

    guide_color = :black
    guide_width = 1
    guide_settings = (color=guide_color, linewidth=guide_width, linestyle=:dash)
    hms = map(zip(1:4, axs)) do (idx, ax)
        _i, _j = to_bin(idx - 1)
        i = 2 - _i
        j = 2 * _j + 1
        hm = heatmap!(
            ax,
            range(1, 4, length=201),
            range(0, 3, length=201),
            (m, c_c) -> nash_equilibrium(m, c_c, 1)[idx],
            colormap=cgrads[idx],
            colorrange=(0, 1),
            # lowclip=:white
        )
        cb = Colorbar(fig[i, j+1], hm, label="Proportion of agents", tellheight=true)
        # cb.height = Relative(0.73)
        vlines!(ax, [2]; guide_settings...)
        ablines!(ax, [-1], [1]; guide_settings...)
        linesegments!(ax, [(Point2(2, 2), Point2(3, 3))]; guide_settings...)
        hm
    end
    Label(fig[0, :], font=:bold, text="Prevalence of each strategy in the Nash Equilibrium", justification=:center)
    for filetype in ("png", "pdf")
        save("figures/nash-equilibrium-mixed-2x2.$filetype", fig)
    end
    fig
end

function u(focal_strategy::Integer, strategies, sp)
    (; c, a, m, ϵD, ϵA) = sp
    # NP, NA, DP, DA
    # [0, 0], [0, 1], [1, 0], [1, 1]
    # 00, 10, 01, 11 (note swapped)
    R_1 = get_payoff_matrix(c, a, m, ϵD, ϵA)
    return sum(R_1[focal_strategy+1, i] * strategies[i] for i in 1:4)
end

