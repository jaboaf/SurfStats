using GRUtils

# Creating a vector field in R^2
V = transpose([repeat(-5:5, inner=11) repeat(-5:5, outer=11)])

I = [1 0; 0 1] # Identity perm
T = [0 1; 1 0] # Transposition perm

IV = mapslices(v-> I*v, V, dims=1)
TV = mapslices(v-> T*v, V, dims=1)

hold(true)
quiver( V[1,:],V[2,:] , IV[1,:],IV[2,:] )
quiver( V[1,:],V[2,:] , TV[1,:],TV[2,:] )
aspectratio(1)
legend("v -> Iv", "v -> Tv")


