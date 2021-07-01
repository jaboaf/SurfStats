# D_n = \sup_x |F_n(x) - F(x)|
include("DataAndDefns.jl")
using GRUtils

B = partitionBy(:m_c,val=:λ_c);
N_B = length.(last.(B))
𝔹 = map(B) do b
		b[1] => 1/length(b[2])*mapreduce(e_G,+,b[2])
	end;
# Bunch of 1s modulo some precision issues.
sum.(last.(𝔹));

nonsym_B = map(𝔹) do 𝕓
	𝕓[1] => abs.(𝕓[2]-𝐒(𝕓[2]))
end;
# a picture of the magnitudes
plot(sum.(last.(nonsym_B)))
# TV should be 1/2 |𝕓 - S(𝕓)|
plot(1/2*sum.(last.(nonsym_B)))
# why dont we sort these to get a better sense of distribution of magnitudes
plot(1/2 * sort(sum.(last.(nonsym_B)),rev=true) )
histogram(1/2 * sum.(last.(nonsym_B)) )
# now see devs as a function of length of process
scatter(N_B,maximum.(last.(nonsym_B)) )
# We can plot upper bound on signifigance
# ie. dots are under curve with probability 1-α
d(a,N) = sqrt(-log(a/(2*factorial(5)))/(2*N))
oplot(LinRange(1,880,880),n -> d(0.01,n), label= "α = 0.01")
oplot(LinRange(1,880,880),n -> d(0.05,n), label= "α = 0.05")
