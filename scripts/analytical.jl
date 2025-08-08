using StaticArrays
using CairoMakie
using JLD2

include("../implementation/payoff_from_interaction.jl")
include("../implementation/common.jl")
include("../implementation/structured.jl")


function u(strategies, parameters)
    (; c_c, c_p, m, ϵ_p, ϵ_c) = parameters
    R_1 = get_payoff_matrix(c_p, c_c, m, ϵ_p, ϵ_c)
    # display(R_1)
    return strategies' * (R_1' * strategies)
end

function average_prevalence(sp, N)
    S_initial = rand_S_initial_structured(sp.Z_1, sp.Z_2)
    strategy_count_by_generation = main_simulation_loop(S_initial, N, sp)
    return [sum(mean.(rows)) for rows in Iterators.partition(eachrow(strategy_count_by_generation), 4)][1:4]
end

function mixed_nash(m, c_c, c_p)
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

function plot_f(new_m, new_a, sp, N)
    new_sp = StructuredParameters(sp.Z_1, sp.Z_2, sp.β, sp.μ, sp.c_p, new_a, new_m, new_m, sp.α, sp.ϵ_p, sp.ϵ_c)
    AP = average_prevalence(new_sp, N) / (sp.Z_1 + sp.Z_2)
    MN = mixed_nash(new_m, new_a, 1)
    Δ = sum(abs, AP .- MN)
    println("($new_m, $new_a): $AP, $MN, $Δ")
    return Δ
end

function plot_u_f(new_m, new_a, sp, N)
    new_sp = StructuredParameters(sp.Z_1, sp.Z_2, sp.β, sp.μ, sp.c_p, new_a, new_m, new_m, sp.α, sp.ϵ_p, sp.ϵ_c)
    AP = average_prevalence(new_sp, N) / (sp.Z_1 + sp.Z_2)
    MN = mixed_nash(new_m, new_a, 1)
    parameters_ntuple = (; c_c=new_a, c_p=sp.c_p, m=new_m, ϵ_p=sp.ϵ_p, ϵ_c=sp.ϵ_c)
    u_AP = u(AP, parameters_ntuple)
    u_MN = u(MN, parameters_ntuple)
    # @show u_AP, u_MN
    # println("$(abs(u_AP - u_MN))")
    # @show u_min
    Δu = abs(u_AP - u_MN)
    println("($new_m, $new_a): $(round.(AP;sigdigits=3)), $(round.(MN;sigdigits=3)), $(round(Δu;sigdigits=3))")
    return Δu
end

function get_AP_MN(new_m, new_a, sp, N)
    println("($new_m, $new_a)")
    new_sp = StructuredParameters(sp.Z_1, sp.Z_2, sp.β, sp.μ, sp.c_p, new_a, new_m, new_m, sp.α, sp.ϵ_p, sp.ϵ_c)
    AP = average_prevalence(new_sp, N) / (sp.Z_1 + sp.Z_2)
    MN = mixed_nash(new_m, new_a, 1)
    return (; AP, MN)
end


begin
    Z_1 = 100
    Z_2 = 0
    β = 1
    μ = 1 / Z_1 + Z_2
    c_p = 1
    c_c = 1
    m_in = 1
    m_out = 1
    α = 1
    ϵ_p = 0.1
    ϵ_c = 0.1
end
sp = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c, m_in, m_out, α, ϵ_p, ϵ_c)
N_diff = 5_000
m_range = range(1, 4, length=101)
c_c_range = range(0, 3, length=101)
data_AP_MN = map(Iterators.product(m_range, c_c_range)) do (m, c)
    get_AP_MN(m, c, sp, N_diff)
end
APs = getindex.(data_AP_MN, :AP)
MNs = getindex.(data_AP_MN, :MN)

# @save "data_AP_MN.jld2" data_AP_MN

get_strategy_diff(AP, MN) = sum(abs(i - j) for (i, j) in zip(AP, MN))
function get_util_diff(AP, MN, (m, c_c), sp)
    parameters_ntuple = (; c_c, m, c_p=sp.c_p, ϵ_p=sp.ϵ_p, ϵ_c=sp.ϵ_c)
    u_AP = u(AP, parameters_ntuple)
    u_MN = u(MN, parameters_ntuple)
    Δu = abs(u_AP - u_MN)
    return Δu
end
get_strategy_diff.(APs, MNs)
get_util_diff.(APs, MNs, Iterators.product(m_range, c_c_range), Ref(sp))


begin
    figsize = (750, 300)
    fig = Figure(; size=figsize)
    ax_stratdiff = Axis(fig[1, 1];
        aspect=1,
        title="Total Absolute Difference in Population State",
        # subtitle="(Simulated vs Analytical)"
    )
    ax_stratdiff.xlabel = "Contribution multiplier"
    ax_stratdiff.ylabel = "Claiming cost"

    ax_utildiff = Axis(fig[1, 3];
        aspect=1,
        title="Absolute Difference in Utility",
        # subtitle="(Simulated vs Analytical)"
    )
    ax_utildiff.xlabel = "Contribution multiplier"
    ax_utildiff.ylabel = "Claiming cost"

    hm_stratdiff = heatmap!(
        ax_stratdiff,
        m_range,
        c_c_range,
        get_strategy_diff.(APs, MNs),
        colormap=:viridis,
        # colorrange=(-1, 1),
        # lowclip=:white,
        # tellheight=true,
    )
    hm_utildiff = heatmap!(
        ax_utildiff,
        m_range,
        c_c_range,
        get_util_diff.(APs, MNs, Iterators.product(m_range, c_c_range), Ref(sp)),
        colormap=:thermal,
        # colorrange=(-1, 1),
        # lowclip=:white,
        # tellheight=true,
    )
    guide_color = :white
    guide_width = 2
    guide_settings = (color=guide_color, linewidth=guide_width, linestyle=:dash)
    for ax in (ax_stratdiff, ax_utildiff)
        vlines!(ax, [2]; guide_settings...)
        ablines!(ax, [-1], [1]; guide_settings...)
        linesegments!(ax, [(Point2(2, 2), Point2(3, 3))]; guide_settings...)
    end
    cb_stratdiff = Colorbar(fig[1, 2], hm_stratdiff, width=10, height=220, tellheight=true, label="Difference in Population State")
    cb_utildiff = Colorbar(fig[1, 4], hm_utildiff, width=10, tellheight=true, label="Difference in Utility")
    # for filetype in ("png", "pdf")
    #     # save("figures/simulation_theoretical_difference.$filetype", fig)
    # end
    display(fig)
end

for filetype in ("png", "pdf")
    # save("figures/simulation_theoretical_difference.$filetype", fig)
end