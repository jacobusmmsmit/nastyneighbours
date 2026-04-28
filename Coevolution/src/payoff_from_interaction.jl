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

function get_payoff_matrix(pot::AbstractVector, c::Real, a::Real)
    #! format: off
    payoff_matrix = SA[
        1/2         0               pot[1]/2    0;
        1-a         1/2-a           pot[1]-a    pot[1]/2 - a;
        pot[1]/2-a  -c              pot[2]/2-c  -c;
        pot[1]-c-a  pot[1]/2-c-a    pot[2]-c    pot[2]/2-c-a
    ]
    #! format on
    return payoff_matrix
    
end

function add_errors_to_payoff_matrix(payoff_matrix, ϵ_p, ϵ_c)
    return @SMatrix [
        begin
            my_errors = SMatrix{4, 4, Float64}(repeat_vector(strategic_errors(i - 1, ϵ_p, ϵ_c)))
            their_errors = SMatrix{4, 4, Float64}(repeat_vector(strategic_errors(j - 1, ϵ_p, ϵ_c)))
            likelihoods = my_errors * (their_errors') / 4
            sum(payoff_matrix .* likelihoods)
        end for i in 1:4, j in 1:4
    ]
end

function payoff_from_interaction(si::SVector{2, Bool}, sj::SVector{2, Bool}, payoff_matrix_with_errors)
    i_claim, i_produce = si
    j_claim, j_produce = sj
    i_idx = 1 + i_claim + 2i_produce
    j_idx = 1 + j_claim + 2j_produce
    return payoff_matrix_with_errors[i_idx, j_idx]
end

get_pot(b0,b1,b2) = SA[b0*b0, b0*(b0+b1), (b0+b2)*(b0+b2)]