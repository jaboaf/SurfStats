using Plots
gr()
using Distributions
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
function wtOfPerm(τ::Array, D::Array)
    p = 1
    for i in 1:length(τ)
        p += D[i]*D[τ[i]]
    end
    return p
end

# Beta given "friendliness" params
# β = [ ... ] means "friendliness" of node i = β[i]
# Technically Not the beta model
function ProbGivenβ(β::Array)
    b = length(β)
    Sb = permGroup(b)
    β /= norm(β,1)
    βresult = map(x->wtOfPerm(x,β), Sb)
    βresult /= sum(βresult)
    return βresult
end

β = [1, 1, 1, 2, 2, 4]
#β .-= mean(β)
β /= norm(β,1)

b = length(β)
U = rand(Normal(1, .5), b)
#U .-= mean(U)
U /= norm(U,1)
Sb = permGroup(b)

Uresult = map(x->wtOfPerm(x,U), Sb)
Uresult /= sum(Uresult)
βresult = map(x->wtOfPerm(x,β), Sb)
βS = sum(βresult)
βresult /= βS
B = β'*β
mapreduce(x->wtOfPerm(x,β)*PMatrix(x),+, Sb) / βS

βresult
histogram(βresult, bins=30)
histogram!(Uresult)
Sb[findall(x-> x>maximum(Uresult), βresult)]

sum(PMatrix(rand(Sb)) .* B)

dadsSample = []
for i in 1:1000
    push!(dadsSample, sum(sample(βresult, 100, replace=false)) )
end
histogram(dadsSample)
β

s = rand(1000)
(s .- mean( s )) / std(s)
histogram(s)

dist0 = ProbGivenβ(rand(Normal(10,1),6))
histogram(dist0)

v = rand(Normal(10,1),7)
dist10 = ProbGivenβ( v )
histogram(dist10)

distNeg10 = ProbGivenβ(-1*v)
histogram!(distNeg10)

histogram()
D = []
for n in 1:8
    push!(D, ProbGivenβ( rand(n) ))
end

histogram(D[2] .- median(D[2]) )
histogram!(D[3] .- median(D[3]) )
histogram!(D[4] .- median(D[4]) )
histogram!(D[5] .- median(D[5]) )
histogram(D[6] .- median(D[6]) )





graphplot([0 1 0; 0 0 1; 1 0 0])
allFours = ProbGivenβ( 4*ones(6) .+ rand(6) )
histogram(allFours)
