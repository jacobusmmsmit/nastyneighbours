using StaticArrays
using Random
using Chairmarks

abstract type AbstractSimulationParameters end

struct StructuredSimulationParameters <: AbstractSimulationParameters
    N::Int
    Z::Int
    Z_1::Int
    Z_2::Int
    β::Float64
    μ::Float64
    c_p::Float64
    c_c::Float64
    m_in::Float64
    m_out::Float64
    α::Float64
    ϵ_p::Float64
    ϵ_c::Float64
end

struct UnstructuredSimulationParameters <: AbstractSimulationParameters
    N::Int
    Z::Int
    β::Float64
    μ::Float64
    c_p::Float64
    c_c::Float64
    m::Float64
    ϵ_p::Float64
    ϵ_c::Float64
end

function main_simulation_loop(S_initial::AbstractMatrix{Bool}, sp::AbstractSimulationParameters)
    (; Z, N, μ, β) = sp
    S = S_initial
    T = zeros(Bool, (size(S, 1), Z, N)) # Strategy by Agents by Generation
    strategy_update_order = collect(1:Z)
    for G in 1:N
        for i in shuffle!(strategy_update_order)
            if rand() < μ
                for k in axes(S, 1)
                    S[k, i] = rand(Bool)
                end
            else
                j = let
                    _j = rand(1:Z-1)
                    _j * (_j < i) + (_j + 1) * (_j ≥ i)
                end
                U_i = average_utility(i, S, sp)
                U_j = average_utility(j, S, sp)
                P_ij = inv(1 + exp(-β * (U_j - U_i)))
                if rand() < P_ij
                    @views S[:, i] .= S[:, j]
                end
            end
        end
        @views T[:, :, G] .= S
    end
    return T
end

function average_utility(i, S, sp::StructuredSimulationParameters)
    (; Z, Z_1, Z_2, α, m_in, m_out, c_p, c_c, ϵ_p, ϵ_c) = sp
    g_i = i > Z_1 # 0 => g1, 1 => g2
    Z_in, Z_out = !g_i .* (Z_1, Z_2) + g_i .* (Z_2, Z_1)
    W = 1 / (α * (Z_in - 1) + (1 - α) * Z_out)
    U_i = 0.0
    for j in 1:Z
        j == i && continue
        g_j = j > Z_1 # 0 => group 1, 1 => group 2
        if g_i == g_j
            w = α
            m = m_in
        else
            w = 1 - α
            m = m_out
        end
        U_i += w / W * payoff_from_interaction(i, j, S, c_p, c_c, m, ϵ_p, ϵ_c)
    end
    return U_i
end

function average_utility(i, S, sp::UnstructuredSimulationParameters)
    (; Z, c_p, c_c, m, ϵ_p, ϵ_c) = sp
    U_i = 0.0
    for j in 1:Z
        j == i && continue
        U_i += payoff_from_interaction(i, j, S, c_p, c_c, m, ϵ_p, ϵ_c) / Z
    end
    return U_i
end

function payoff_from_interaction(i, j, S, c_p, c_c, m, ϵ_p, ϵ_c)
    @views i_produce, i_claim = S[1:2, i] # adjust for which group (1:2, or 3:4)
    @views j_produce, j_claim = S[1:2, j] # adjust for group
    common_resource = (i_produce + j_produce) * c_p * m
    i_partition = (i_claim & !j_claim) + 0.5(!⊻(i_claim, j_claim))
    return common_resource * i_partition - i_claim * c_c
end


N = 10_000
Z = 50
β = 1
μ = 1/Z
c_p = 1
c_c = 1.5
m = 2.5
ϵ_p = 0
ϵ_c = 0

usp = UnstructuredSimulationParameters(N, Z, β, μ, c_p, c_c, m, ϵ_p, ϵ_c)
S_initial = rand(Bool, 2, Z)

@btime main_simulation_loop($S_initial, $usp)