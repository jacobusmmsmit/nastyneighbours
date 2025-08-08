include("../implementation/structured.jl")
include("../implementation/unstructured.jl")

Z_1 = 50
Z_2 = 0
Z = Z_1 + Z_2
S_structured = rand_S_initial_structured(Z_1, Z_2)
S_unstructured = let
    _S = [sum(rows) for rows in Iterators.partition(S_structured, 4)]
    _S[1:4]
end

begin
    β = 1
    μ = 1 / (Z_1+Z_2)
    c_p = 1
    c_c = 1
    m_in = 1.5
    m_out = 1.5
    m = 1.5
    α = 0.5
    ϵ_p = 0
    ϵ_c = 0
    N = 10_000
end

sp = StructuredParameters(Z_1, Z_2, β, μ, c_p, c_c, m, m, α, ϵ_p, ϵ_c)
up = UnstructuredParameters(Z_1 + Z_2, β, μ, c_p, c_c, m, ϵ_p, ϵ_c)

structured_strategies = vcat(strategies, strategies)
for _ in 1:10
    println("------------------------------")
    gi = true
    fw = FrequencyWeights(S_structured)
    si_structured_idx = sample(1:32, fw)
    si_structured = structured_strategies[si_structured_idx]
    # @show si_structured_idx, si_structured, S_structured[si_structured_idx]
    si_unstructured = SVector{2,Bool}(si_structured[3], si_structured[4])
    aus = average_utility(si_structured, false, S_structured, sp)
    auu = average_utility(si_unstructured, S_unstructured, up)
    if !(aus≈ auu)
        println("($aus, $auu)")
    end
end