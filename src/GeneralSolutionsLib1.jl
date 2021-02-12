using Plots
gr()
using Distributions
using Statistics
using LinearAlgebra

# Distribution on permutations given a degree sequence
function permGroup(A::Union{Set,Array})
    unique!(A)
    P = Array{eltype(A),1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(P, head)
        end
    end
    continuePerm(eltype(A)[], A)
    return P
end
function permGroup(n::Integer)
    if n > 14 error(" you gave $n .... thats $(factorial(n)) element") end
    A = collect(Int8, 1:n)
    P = Array{eltype(A),1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(P, head)
        end
    end
    continuePerm(eltype(A)[], A)
    return P
end


# Suppose T(g) = degree sequence of g
# t = sort(rand(1:5,5), rev=true)
t = [2,1,1,2,4,1,2,1]
function graphFromDeg(T::Array)
    G = Array[]
    t = sort(copy(T), rev=true)
    dotDiagram(t)
    println("This is the degree sequence you passed (yes, we sorted it)")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    while any(t .> 0)
        rem = findall(n-> n>0, t)
        m = minimum(t[rem])
        gs = rand(permGroup(rem), m)
        append!(G, gs)
        dotDiagram(t, forceMin=1)
        println("------------------------")
        println("sampled $m perms from Sym($rem)")
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        t .-= m
    end
    dotDiagram(t, forceMin=1)
    println("------------------------")
    println("Nothing left? Correct.")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    println("Now lets see whats under the hood: \n")
    dotDiagram(t)
    return G
end
#=
More succinct version of graphFromDeg

function graphFromDeg(t::Array)
    G = []
    t = copy(T)
    while any(t .> 0)
        rem = findall(n-> n>0, t)
        m = minimum(t[rem])
        gs = rand(permGroup(rem), m)
        append!(G, gs)
        global t[rem] .-= m
    end
    return G
end

=#

tt = [1,3,2,2,1,4,3,2]
A = zeros(Int8, (length(tt),maximum(tt)))
for (i,v) in enumerate(tt)
    A[i,1:v] .= 1
end

sort!(tt, rev=true)
B = zeros(Int8, (length(tt),maximum(tt)))
for (i,v) in enumerate(tt)
    B[i,1:v] .= 1
end



function dotDiagram(A::Array; forceMin=nothing)
    if forceMin != nothing min=forceMin
    else min = minimum(A)
    end
    for d in maximum(A):-1:min
        print(d, "\t")
        for a in A
            if a>=d print("•  ")
            else print("   ")
            end
        end
        print("\n")
    end
    println("V","\t", join(1:length(A), "  "))
end

function makeAdj(H::Array{Array,1})
    n = maximum(length.(H))
    A = zeros(Int8, n,n)
    for τ in H
        k = length(τ)
        for i in 1:k
            A[i,τ[i]] = 1
        end
        #A -= diagm([ones(Int8,k)... , zeros(Int8, n-k)...])
    end
    return A
end

function wtPMatrix(τ::Array, D::Array)
    p = zeros(n,n)
    for i in 1:n
        p[ i , τ[i] ] = D[i]*D[τ[i]]
    end
    return p
end

function PMatrix(τ::Array; inSnWithn=nothing)
    if inSnWithn==nothing
        p = zeros(Int16, length(τ),length(τ))
    else
        p = zeros(Int16, inSnWithn,inSnWithn)
    end
    for i in 1:length(τ)
        p[ i , τ[i] ] = 1
    end
    return p
end


# Beta given "friendliness" params
# β = [ ... ] means "friendliness" of node i is β[i]
# THIS IS NOT THE ACTUAL BETA MODEL
# actual beta model lets
function ProbGivenβ(β::Array; rep=:perm)
    n = length(β)
    Sn = permGroup(n)
    function wtOfPerm(τ::Array, D::Array)
        p = 1
        for i in 1:n
            p += D[i]*D[τ[i]]
        end
        return p
    end
    function wtPermMatrix(τ::Array, D::Array)
        p = zeros(n,n)
        for i in 1:n
            p[ i , τ[i] ] = D[i]*D[τ[i]]
        end
        return p
    end
    if rep === :perm
        βresult = map(x->wtOfPerm(x,β), Sn)
        βresult /= sum(βresult)
    elseif rep === :matrix
        βresult = sum(map(x->wtPermMatrix(x,β), Sn))
        βresult *= sum(βresult)^-1
    else
        error(" bad arg for rep, pass :perm or :matrix")
    end
    return βresult
end

β = [rand(Normal(4,1),5)... ,6]
b = length(β)
U = mean(β)*ones(b) + 1*rand(b)
Sb = permGroup( collect( 1:b ))
function wtOfPerm(τ::Array, D::Array)
    p = 1
    for i in 1:b
        p += D[i]*D[τ[i]]
    end
    return p
end
Uresult = map(x->wtOfPerm(x,U), Sb)
Uresult /= sum(Uresult)
βresult = map(x->wtOfPerm(x,β), Sb)
βresult /= sum(βresult)

histogram(βresult, bins=50)
histogram!(Uresult)
Sb[findall(x-> x>maximum(Uresult), βresult)]

dadsSample = []
for i in 1:1000
    push!(dadsSample, sum(sample(βresult, 100, replace=false)) )
end
histogram(dadsSample)
β

s = rand(1000)
(s .- mean( s )) / std(s)
histogram(s)

dist0 = ProbGivenβ(rand(Normal(0,1),6))
histogram(dist0)

dist10 = ProbGivenβ(rand(Normal(10,1),6))
histogram!(dist10)

distNeg10 = ProbGivenβ(rand(Normal(-10,1),6))
histogram!(distNeg10)







graphplot([0 1 0; 0 0 1; 1 0 0])
allFours = ProbGivenβ( 4*ones(6) .+ rand(6) )
histogram(allFours)
