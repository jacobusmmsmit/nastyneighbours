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


# function payoff_from_interaction_old(si::SVector{2,Bool}, sj::SVector{2,Bool}, c_p::Number, c_c::Number, m::Number, ϵ_p::Number, ϵ_c::Number)
#     i_claim, i_produce = si
#     j_claim, j_produce = sj
#     common_resource = (i_produce + j_produce) * c_p * m
#     i_partition = (i_claim & !j_claim) + 0.5(!xor(i_claim, j_claim))
#     return common_resource * i_partition - (i_produce * c_p) - (i_claim * c_c)
# end

function get_payoff_matrix(pot::AbstractVector, c::Float64, vs::AbstractVector{Float64}, as::AbstractVector{Float64})
    vN, vC, vP, vPC = vs
    aN, aC, aP, aPC = as
    #! format: off
    payoff_matrix = SA[
        pot[1]/2                     1/2*(1 - vN)*pot[1]                 pot[2]/2                     1/2*(1 - vN)*pot[2];
        -aN + 1/2*(1 + vN)*pot[1]    -aC + pot[1]/2                      -aP + 1/2*(1 + vP)*pot[2]    -aPC + 1/2*(1 - vC + vPC)*pot[2];
        -c + pot[2]/2                -c + 1/2*(1 - vP)*pot[2]            -c + pot[3]/2                -c + 1/2*(1 - vP)*pot[3];
        -aN - c + 1/2*(1 + vN)*pot[2]  -aC - c + 1/2*(1 + vC - vPC)*pot[2]  -aP - c + 1/2*(1 + vP)*pot[3]  -aPC - c + pot[3]/2
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