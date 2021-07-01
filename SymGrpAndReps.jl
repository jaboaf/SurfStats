function Sym(A::Array)
    if length(A) >= 14
        error("Slow down... be careful, n >= 14, and |S_n| = $(factorial(big(length(A)))) ")
    end
    Perms = Array{eltype(A),1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for (i,t) in enumerate(tail)
                newHead = vcat(head, t)
                newTail = vcat(tail[1:(i-1)], tail[(i+1):end])
                continuePerm(newHead, newTail)
            end
        else
            push!(Perms, head)
        end
    end
    continuePerm(eltype(A)[], A)
    return Perms
end

function Sym(A::Set)
    if length(A) >= 14
        error("Slow down... be careful, n >= 14, and |S_n| = $(factorial(big(length(A)))) ")
    end
    Perms = Array{eltype(A),1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(Perms, head)
        end
    end
    continuePerm(eltype(A)[], A)
    return Perms
end

function Sym(n::Integer)
    if n >= 14
        error("Slow down... be careful, n >= 14, and |S_n| = $(factorial(big(n))) ")
    end
    A = collect(Int8, 1:n)
    Perms = Array{Int8,1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(Perms, head)
        end
    end
    continuePerm(Int8[], A)
    return Perms
end

function RepSym(n::Integer)
    if n >= 14
        error("Slow down... be careful, n >= 14, and |S_n| = $(factorial(big(n))) ")
    end
    MatPerms = BitArray{2}[]
    function continuePerm( M ,tail)
        if length(tail) > 0
            c = n - length(tail) + 1
            for t in tail
                newM = deepcopy(M)
                newM[ t , c ] = 1
                newTail = setdiff(tail, [t])
                continuePerm( newM , newTail)
            end
        else
            push!( MatPerms, M )
        end
    end
    # Init function w/ nxn BitArray of 0s and dense  [1,2,...,n]
    continuePerm( falses(n,n) , BitSet(1:n) )
    return MatPerms
end

function Rep(t::Array)
    M = falses(length(t),length(t))
    for i in 1:length(t)
        M[ t[i] , i ] = 1
    end
    return M
end

function Rep(t::Array, Res::Array)
    if !(issubset(Res, t))
        error("Restriction is not contained in the domain of the perm, t.")
    end
    M = falses(length(t), length(t))
    for j in Res
        M[ t[j] , j ] = 1
    end
    return M
end

function Rep(t::Array, Ind::Array)
    if !(issubset(t, Ind))
        error("The image of the perm t is not contained in the Induced space")
    end
    M = falses( length(Ind), length(Ind))
    for j in sort(t)
        M[ t[j] , j ] = 1
    end
    return M
end

#function TensorRep(t::)


function sgn(a::Array{Int8,1})
    d = length(a)
    return prod([ sign(a[j] - a[i]) for i in 1:d for j in i+1:d ])
end

function sgn(M::BitArray{2})
    d = size(M,1)
    negs = 0
    for i in 1:d
        ti = findfirst( M[:,i] )
        negs += sum( M[ 1:(ti-1) , (i+1):d ] )
    end
    return negs%2 == 0 ? 1 : -1
end


# Symmetrization operations
# for a permutation
function 𝐒( M::Array{<:Number, d}) where d
    return 1//factorial(d) * mapreduce( p -> permutedims(M, p), +, Sym(d) )
end
# for a d-dimensional array
function 𝐒( M::BitArray{d}) where d
    return 1/factorial(d) * mapreduce( p -> permutedims(M, p), +, Sym(d) )
end

# Alternating operations
# for a permutation
function 𝐀( M::Array{<:Number,d}) where d
    return 1//factorial(d)* mapreduce( p -> sgn(p)*permutedims(M, p), +, Sym(d) )
end
# for a d-dimensional array
function 𝐀( M::BitArray{d}) where d
    return 1/factorial(d)* mapreduce( p -> sgn(p)*permutedims(M, p), +, Sym(d) )
end

#=
setindex!(::typeof(f), ::Int64, ::Int64)

f(1) = 1
struct Group <: Union{ Set{Function}, Function }
A::Set{Function}
op::Function
assert( all( g * h for a in A, h in A))
end

struct Group <: Union{ Set{Function}, Function }
A::Set{Function}
op::Function
assert( all( g * h for a in A, h in A))
end

Type Group <: Union{Set,Function} end

struct Group


abstract type Group <: Union{ Set{T}, Function{T,T}{T} } end


Sym(C::Set) = ({ All Permuations of C } , ∘ )
Group = ( {things} , op: )

struct Group{T}
    S::Set{T}
    op::Function
end

struct SymmetricGroup
    C::Set{Countries}
    op::Function{Countries,Countries}
end
    
=#







