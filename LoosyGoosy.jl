include("DataAndDefns.jl")
L = map(x->x.λ_c,WAVES);
K = vcat(collect.(L)...);
