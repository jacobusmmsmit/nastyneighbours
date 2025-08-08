one_at(i, n) = (x = zeros(n); x[i] = 1.0; return x)
to_bin(i) = SVector{2,Bool}(i >> shift & 1 != 0 for shift in false:true)

function strategic_errors(s::SVector{2,Bool}, ϵD, ϵA)
    bins = SVector{4}(to_bin(i) for i in 0:3)
    # NP, NA, DP, DA
    return map(bins) do b
        x = first(b) == first(s) ? 1 - ϵA : ϵA
        y = last(b) == last(s) ? 1 - ϵD : ϵD
        x * y
    end
end

function strategic_errors(strat::Integer, ϵD, ϵA)
    return strategic_errors(to_bin(strat), ϵD, ϵA)
end

function repeat_vector(v::SVector{N,T}) where {N,T}
    vov = SVector{N,SVector{N,T}}(v for _ in 1:N)
    return reinterpret(reshape, T, vov)
end


function payoff_from_interaction_old(si::SVector{2, Bool}, sj::SVector{2, Bool}, c_p::Number, c_c::Number, m::Number, ϵ_p::Number, ϵ_c::Number)
    i_claim, i_produce = si
    j_claim, j_produce = sj
    common_resource = (i_produce + j_produce) * c_p * m
    i_partition = (i_claim & !j_claim) + 0.5(!xor(i_claim, j_claim))
    return common_resource * i_partition - (i_produce * c_p) - (i_claim * c_c)
end

function get_payoff_matrix(c_p, c_c, m, ϵ_p, ϵ_c)
    # full_payoff_matrix = SA[
    #     0          0            m*c/2       0;
    #     -a         -a           m*c-a        m*c/2-a;
    #     (m*c)/2-c -c           m*(2c)/2-c -c;
    #     m*c-c-a    (m*c)/2-c-a m*(2c)-c-a  m*(2c)/2-c-a
    # ]
    R = SA[
        0 0 m*c_p/2 0;
        -c_c -c_c m*c_p-c_c m*c_p/2-c_c;
        (m*c_p)/2-c_p -c_p m*(2c_p)/2-c_p -c_p;
        m*c_p-c_p-c_c (m*c_p)/2-c_p-c_c m*(2c_p)-c_p-c_c m*(2c_p)/2-c_p-c_c
    ]
    R̃ = @SMatrix [
        begin
            my_errors = SMatrix{4,4,Float64}(repeat_vector(strategic_errors(i - 1, ϵ_p, ϵ_c)))
            their_errors = SMatrix{4,4,Float64}(repeat_vector(strategic_errors(j - 1, ϵ_p, ϵ_c)))
            likelihoods = my_errors * (their_errors') / 4
            sum(R .* likelihoods)
        end for i in 1:4, j in 1:4
    ]
    return R̃
end

function payoff_from_interaction(si::SVector{2, Bool}, sj::SVector{2, Bool}, c_p::Number, c_c::Number, m::Number, ϵ_p::Number, ϵ_c::Number)
    payoff_matrix = get_payoff_matrix(c_p, c_c, m, ϵ_p, ϵ_c)
    i_claim, i_produce = si
    j_claim, j_produce = sj
    i_idx = 1 + i_claim + 2i_produce
    j_idx = 1 + j_claim + 2j_produce
    return payoff_matrix[i_idx, j_idx]
end