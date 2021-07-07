include("DataAndDefns.jl")
L = map(x->x.λ_c,WAVES);
K = vcat(collect.(L)...);

𝕐 = e_G.(L);
𝐒𝕐 = 𝐒.(𝕐);
𝕖 = 𝕐 - 𝐒𝕐