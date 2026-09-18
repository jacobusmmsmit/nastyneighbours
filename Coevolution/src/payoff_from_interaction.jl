function strategic_errors(s::SVector{2,Bool}, ϵD, ϵA)
    function inner(b)
        x = first(b) == first(s) ? 1 - ϵA : ϵA
        y = last(b) == last(s) ? 1 - ϵD : ϵD
        return x * y
    end
    bins = SVector{4,SVector{2,Bool}}(SA[0, 0],
        SA[1, 0],
        SA[0, 1],
        SA[1, 1]
    )
    # NP, NA, DP, DA
    return inner.(bins)
end

function strategic_errors(strat::Integer, ϵD, ϵA)
    return strategic_errors(to_bin(strat), ϵD, ϵA)
end

function strategic_error_matrix(ϵD, ϵA)
    SMatrix{4,4,Float64}(
        # strategy 00
        (1-ϵA) * (1-ϵD), ϵA*(1-ϵD), (1-ϵA)*ϵD, ϵA*ϵD,

        # strategy 10
        ϵA * (1-ϵD), (1-ϵA)*(1-ϵD), ϵA*ϵD, (1-ϵA)*ϵD,

        # strategy 01
        (1-ϵA) * ϵD, ϵA*ϵD, (1-ϵA)*(1-ϵD), ϵA*(1-ϵD),

        # strategy 11
        ϵA * ϵD, (1-ϵA)*ϵD, ϵA*(1-ϵD), (1-ϵA)*(1-ϵD)
    )
end

function get_error_likelihoods(ϵD, ϵA)
    errors = strategic_error_matrix(ϵD, ϵA)
    return errors * errors' / 4
end

function get_payoff_matrix(pot::AbstractVector, c::Real, a::Real, aggression_ratio::Rational=1 // 2)
    #! format: off
    payoff_matrix = SA[
        1/2         0               pot[1]/2    0;
        1-a         aggression_ratio-a           pot[1]-a    pot[1]*aggression_ratio - a;
        pot[1]/2-c  -c              pot[2]/2-c  -c;
        pot[1]-c-a  pot[1]*aggression_ratio-c-a    pot[2]-c    pot[2]*aggression_ratio-c-a
    ]
    #! format on
    return payoff_matrix
end

function add_errors_to_payoff_matrix(payoff_matrix, ϵ_p, ϵ_c)
    return SMatrix{4, 4, Float64, 16}(
        begin
            my_errors = SMatrix{4, 4, Float64}(repeat_vector(strategic_errors(i - 1, ϵ_p, ϵ_c)))
            their_errors = SMatrix{4, 4, Float64}(repeat_vector(strategic_errors(j - 1, ϵ_p, ϵ_c)))
            likelihoods = my_errors * (their_errors') / 4
            sum(payoff_matrix .* likelihoods)
        end for i in 1:4, j in 1:4
    )
end

function add_errors_to_payoff_matrix(payoff_matrix, likelihoods)
    return SMatrix{4,4,Float64,16}(sum(payoff_matrix .* likelihoods[i,j]) for i in 1:4, j in 1:4)
end

function payoff_from_interaction(si::SVector{2, Bool}, sj::SVector{2, Bool}, payoff_matrix_with_errors)
    i_claim, i_produce = si
    j_claim, j_produce = sj
    i_idx = 1 + i_claim + 2i_produce
    j_idx = 1 + j_claim + 2j_produce
    return payoff_matrix_with_errors[i_idx, j_idx]
end

# get_pot(b0,b1,b2) = SA[b0*b0, b0*(b0+b1), (b0+b2)*(b0+b2)]