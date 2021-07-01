# D_n = \sup_x |F_n(x) - F(x)|
include("DataAndDefns.jl")
using GRUtils

# RV is X = ( λ_c = [ x.λ_c for x in WAVES ])
Y= [w.λ_c for w in WAVES]
n= length(Y)

# We would write..
# 𝔽 = 1/n * ∑ e_{G yᵢ} = 1/n * ∑ 1_{G yᵢ}
𝔽 = 1/n * mapreduce(e_G, +, Y);

F = 𝐒(𝔽);
# The random variable is λ_c

D = abs.(𝔽 - F);
KS = maximum(D)
α = 0.01
sqrt(-log(α/(2*count(!=(0),D)))/(2*n))
sqrt(-log(α/(2*7^5) )/(2*n))
count(!=(0),D)*(n+1)*exp(-2*n*α^2)
7^5 *(n+1)*exp(-2*n*α^2)

d(a) = sqrt(-log(2*a/count(!=(0),D))/(2*n))
d(a,N) = sqrt(-log(a/(2*count(!=(0),D)))/(2*N))

plot(vec(D))
plot(sort(vec(D)))
U = collect( (1/7^5):(100/7^5):1)
N = collect(1:15:n)

plot(U,d.(U))
plot(N,d.(α,N))
surface(N,U,(x,y)->d(y,x))

N_fifty = collect(1:5)
wireframe(N_fifty,U,(x,y)->d(y,x))

