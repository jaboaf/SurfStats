include("SymGrpAndReps.jl")
using GRUtils

# Creating a vector field in R^2
V = transpose([repeat(-5:5, inner=11) repeat(-5:5, outer=11)])

I = [1 0; 0 1] # Identity perm
T = [0 1; 1 0] # Transposition perm

IV = mapslices(v-> I*v, V, dims=1)
TV = mapslices(v-> T*v, V, dims=1)

subplot(1,2, (1))
quiver( V[1,:],V[2,:] , IV[1,:],IV[2,:] , "o" , markersize=.2 )
title("Identity Permuation")
subplot(1,2, (2))
quiver( V[1,:],V[2,:] , TV[1,:],TV[2,:] , "o" , markersize=.2 )
title("Transposition")


# Making one in R^3
Perms = Rep.( Sym(3) )

Pts = Array{Float64}(undef,3, 7*7*7)
x = 1
for i in -3:3; for j in -3:3; for k in -3:3
	Pts[:,x] .= [i;j;k]
	x+=1
end end end

Fields = []
for T in Perms
	push!(Fields, mapslices(z -> T*z, Pts, dims=1) )
end

for i in 1:6
	subplot(3,2, (i) )
	quiver3( Pts[1,:],Pts[2,:],Pts[3,:] , Fields[i][1,:],Fields[i][2,:],Fields[i][3,:], "o", arrowsize=.8, markersize=.2)
end
