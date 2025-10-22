one_at(i, n) = (x = zeros(n); x[i] = 1.0; return x)

to_bin(i) = SVector{2, Bool}(i >> shift & 1 != 0 for shift in false:true)

function repeat_vector(v::SVector{N, T}) where {N, T}
    vov = SVector{N, SVector{N, T}}(v for _ in 1:N)
    return reinterpret(reshape, T, vov)
end
