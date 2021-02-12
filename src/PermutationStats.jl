using Plots; gr()
using Distributions
using Statistics
using LinearAlgebra
using JSON
using Iterators

# This is a helper function that will be useful going forward...
# It constructs SymmetricGroup( A ), and the elements of A can be anything
# note: S_k = permGroup( collect(1:k) )
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

# this constructs a matrix for a given permutation
# inSnWithn
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

# helper function for checking if it is valid to say X < Y
# note: differs from julia's built in issorted which seems to use partial orders
#   ex: issorted([1,2,2,3]) returns true
function isordered(X::Array{T},Y::Array{T}) where T <: Number
    if all(map(t-> isless(t[1],t[2]), Base.product(X,Y)))
        return true
    else
        return false
    end
end

# this tests if an array is strictly ordered by <, which is defined for various things
function strictOrder(S)
    for i in 1:(length(S)-1)
        if !(S[i] < S[i+1])
            return false
        end
    end
    return true
end

# This makes dealing with data quite a bit easier,
#   though the use of symbols can cause difficulties & headaches
isoDict = Dict([  "Australia" => :AUS,
                "Basque Country" => :ESP,
                "Brazil" => :BRA,
                "Fiji" => :FJI,
                "France" => :FRA,
                "Hawaii" => :USA,
                "Indonesia" => :IDN,
                "Italy" => :ITA,
                "Japan" => :JPN,
                "New Zealand" => :NZL,
                "Portugal" => :PRT,
                "South Africa" => :ZAF,
                "Spain" => :ESP,
                "United States" => :USA ])

isoS = collect(values(isoDict)) # For convenience
e = sort(unique(isoS)) # This is our identity permutation

data = Dict()
# File Path for Data
open("Data/CombinedCountries/CleanAllDataCC.txt", "r") do f
    global data
    data = JSON.parse(f)  # parse and transform data
end

waves = []
for wid in keys(data)
    if data[wid]["nJudOrigs"] == 5 & data[wid]["nSubScos"] == 5
        origs = unique(data[wid]["subScoOrig"])
        matchIndicator = (data[wid]["athOrig"] in origs)
        labeledScos = Dict([isoDict[origin] => Float16[] for origin in origs])
        origScoPairs = collect(zip(data[wid]["subScoOrig"],data[wid]["subSco"]))

        labeledScosBinary = Dict([:Match => Float16[], :NoMatch => Float16[] ])

        for p in origScoPairs
            # push!( array of judge scores from country p[1], score=p[2] )
            push!(labeledScos[ isoDict[p[1]] ], p[2])
            if p[1] == data[wid]["athOrig"]
                push!(labeledScosBinary[:Match], p[2])
            else
                push!(labeledScosBinary[:NoMatch], p[2])
            end
        end

        x = (   id=wid,
                evtYear=data[wid]["evtYear"],
                evtOrig=isoDict[data[wid]["evtOrig"]],
                evtName=data[wid]["evtName"],
                evtId=data[wid]["evtId"],
                rnd=data[wid]["rnd"],
                rndId=data[wid]["rndId"],
                heat=data[wid]["heat"],
                heatId=data[wid]["heatId"],
                athName=data[wid]["athName"],
                athId=data[wid]["athId"],
                athOrig=isoDict[data[wid]["athOrig"]],
                currentPoints=data[wid]["currentPoints"],
                endingPoints=data[wid]["endingPoints"],
                panel=labeledScos,
                panelBinary=labeledScosBinary,
                subScos=data[wid]["subSco"],
                subScoOrigs=map(x->isoDict[x], data[wid]["subScoOrig"]),
                panelOrigs=Set(map(x->isoDict[x], data[wid]["subScoOrig"])),
                match=matchIndicator )

        push!(waves, x)
    end
end

BinaryPanels = map(x-> x.panelBinary, filter(w-> all(length.(values(w.panelBinary)).>0), waves) )
length(BinaryPanels)
labels = ["Match", "NoMatch"]
noOrder = 0
Wbinary = []
D = zeros(2,2)
for bp in BinaryPanels
    if isordered(bp[:Match],bp[:NoMatch])
        D += [1 0; 0 1]
        push!(Wbinary, [1 0; 0 1])
    elseif isordered(bp[:NoMatch], bp[:Match])
        D += [0 1; 1 0]
        push!(Wbinary, [0 1; 1 0])
    else
        noOrder += 1
    end
end
# What our result looks like
D
heatmap(labels,labels, D)

# Normalized D
normalizedD = D/length(Wbinary)
heatmap(labels,labels, D)

# Expected Value given uniform prior
unifBinary = [ 1/2 1/2 ; 1/2 1/2]
expected = unifBinary ^ length(Wbinary)


n = length(Wbinary)
s2 = permGroup([1,2])
result = []
for s in 1:500
    Samp = n^-1 * sum(PMatrix.(rand(s2, n)))
    push!(result, maximum(Samp - [1/2 1/2; 1/2 1/2] ) )
end
histogram(result)

sinh(normalizedD)
cosh(normalizedD)


#-----
C = [zeros(2,2) for i in 1:2]
E = zeros(2,2)
for n in 1:20
    E += (length(Wbin))^(-1) * sum( (-1*Wbin) .^ n) / factorial(n)
    C[ (n%2)+1 ] += (length(Wbin))^(-1) * sum( (Wbin) .^ n) / factorial(n)
end
E
C[2] - C[1]^2


for country in isoS
    BinaryPanelByAthOrig = map(w-> w.panelBinary, filter(x-> x.athOrig == country, stdBinaryPanels))
    noOrderBinary = 0
    yesOrderBianry = 0
    DBinary = zeros(2,2)
    for pb in BinaryPanelByAthOrig
        if isordered(pb[:Match],pb[:NoMatch])
            DBinary += [1 0; 0 1]
            yesOrderBianry += 1
        elseif isordered(pb[:NoMatch], pb[:Match])
            DBinary += [0 1; 1 0]
            yesOrderBianry += 1
        else
            noOrderBinary += 1
        end
    end
    println(country)
    println(DBinary)
    println("$yesOrderBianry waves were consistent with a total order {:Match, :NoMatch}")
    println("$noOrderBinary waves were not consistent with any total order on {:Match, :NoMatch}")
    println("----------------------")
end

allPanelCompositions = unique( map(x->x.panelOrigs, waves))
judOrigs = sort(collect(∪(allPanelCompositions...)))
permGroup(judOrigs)
densityOp = Dict([p=>0 for p in Base.product(judOrigs,judOrigs)])
NatToIndex = Dict([orig => i for (i,orig) in enumerate(judOrigs) ])
W = []
consWave = []
totalConsistentPanels = 0
for panelComp in allPanelCompositions
    consistentPanels = 0
    notOnPanel = map(x->NatToIndex[x], setdiff(judOrigs,panelComp) )
    #unifM = zeros(7,7)
    #unifM[notOnPanel,notOnPanel] .= 1 / length(notOnPanel)
    for w in filter(x -> x.panelOrigs == panelComp, waves)
        e = sort(collect(panelComp))
        waveMatrix = zeros(Int,7,7)
        for ord in permGroup(collect(panelComp))
            S = Iterators.product([w.panel[c] for c in ord]...)
            if all(map(strictOrder, S) )
                for (i,ei) in enumerate(e)
                    waveMatrix[ NatToIndex[ei], NatToIndex[ ord[i]] ] = 1
                end
                #waveMatrix += unifM
                consistentPanels +=1
                push!(consWave, w)
                push!(W, waveMatrix)
            end
        end
    end
    totalConsistentPanels += consistentPanels
    println(panelComp)
    println(consistentPanels)
end
println(totalConsistentPanels)
W

maximum(tr.(W))
S = PMatrix.(rand( permGroup(collect(1:7)), length(W)))

Moments = [ (length(W))^-1 * reduce(+, (W) .^ k ) for k in 1:13 ]
heatmap(string.(judOrigs), string.(judOrigs), Moments[4] )

eW = (length(W))^-1 * mapreduce(x-> exp(-1)exp(x)-diagm(ones(7)),+,W)
sum(eW)
heatmap(eW)

eM1 = exp(M1) - diagm(ones(7))
heatmap(eM1)

chi = sum(tanh.(W)) / length(W)
heatmap(chi)

sum(T)
sum([T[i,i] for i in 1:7])
heatmap(T)

LinIndSet = Set(deepcopy(W))
mapreduce(x->x+x',+,LinIndSet)

N = 13 # = ℯ^(7/ℯ)... the upper bound on
D = []
P = zeros(13,13)
for d in 1:13
    C = [zeros(Float16,7,7) for a in 1:d ]
    for n in 1:20
        k = rem(n,d)
        C[ k+1 ] += sum( W .^ (n)) / factorial(n)
    end
    push!(D,C)
    P[1:d,d] .= vec(sum.(C))
end

for i in 1:13
    println("d = $i which gives $(sum.(D[i])) ")
end

heatmap(P, xlabel=" values of d ", ylabel= " values of k (note: k <= d)")
plot(P, xlabel= "n", ylabel= "Total Variation")
for i in 1:13 P[:,i] /= sum(P[:,i]) end

plot(P[:,1], xlabel= "n", xticks=0:1:13, ylabel= "Total Variation")
plot(P[:,1:2], xlabel= "n", xticks=0:1:13, ylabel= "Total Variation")
plot(P[:,1:6], xlabel= "n", xticks=0:1:13, ylabel= "Total Variation")

S = []
Q = zeros(13,13)
for d in 1:13
    C = zeros(Float16,d)
    for n in 1:20
        k = rem(n, d )
        C[ k+1 ] += tr(sum(W .^n)) / factorial(n)
    end
    push!(S,C)
    Q[1:d,d] .= vec(C)
end

length(W)^-1 * sum(tanh.(W))
r = rand(permGroup(collect(1:7)))
tanh(1/7*ones(7,7))


heatmap(string.(collect(1:13)), string.(collect(1:13)),Q, xlabel=" values of d ", ylabel= " values of k (note: k <= d)")
plot(Q, xlabel= "n", ylabel= "Total Variation")
for i in 1:13 Q[:,i] /= sum(abs.(Q[:,i])) end
plot(Q[:,1], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="1st order")
plot!(Q[:,2], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="2nd order")
plot!(Q[:,3], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="3rd order")
plot!(Q[:,4], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="4th order")
plot!(Q[:,5], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="5th order")
plot!(Q[:,6], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="6th order")

plot!([1/factorial(n) for n in 0:13], label="benchmark is e^(n/e)!")


n = 2
plot(x->x, x->1/sum([x^k * Q[n,k] for k in 1:n]), Int(1),Int(13), ylims=(0,1) )
plot!(x->x, x->1/sum([x^k * Q[3,k] for k in 1:3]), Int(1),Int(13), ylims=(0,1) )
plot!(x->x, x->1/sum([x^k * Q[4,k] for k in 1:4]), Int(1),Int(13), ylims=(0,1) )


# ------ Now with fake data ----------

fakeS = []
fakeWaves = [PMatrix.(rand(permGroup(collect(1:7)), 550)); [PMatrix([2,3,1,4,5,6,7]) for i in 1:550] ]
fakeQ = zeros(13,13)
for d in 1:13
    C = zeros(Float16,d)
    for n in 1:20
        k = rem(n, d )
        C[ k+1 ] += tr(sum( (-1)^n *fakeWaves .^n)) / factorial(n)
    end
    push!(fakeS,C)
    fakeQ[1:d,d] .= vec(C)
end
heatmap(string.(collect(1:13)), string.(collect(1:13)),fakeQ, xlabel=" values of d ", ylabel= " values of k (note: k <= d)")
plot(fakeQ, xlabel= "n", xticks=1:13, ylabel= "Total Variation")
for i in 1:13 Q[:,i] /= sum(abs.(Q[:,i])) end
plot(fakeQ[:,1], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="1st order")
plot!(fakeQ[:,2], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="2nd order")
plot!(fakeQ[:,3], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="3rd order")
plot!(fakeQ[:,4], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="4th order")
plot!(fakeQ[:,5], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="5th order")
plot!(fakeQ[:,6], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="6th order")

# ------ Now with symmetric Group on 1... 7 ----------
symS = []
symW = PMatrix.(permGroup(collect(1:7)))
symQ = zeros(13,13)
for d in 1:13
    C = zeros(Float16,d)
    for n in 1:20
        k = rem(n, d )
        C[ k+1 ] += tr(sum( (-1)^n * symW .^n)) / factorial(n)
    end
    push!(symS,C)
    symQ[1:d,d] .= vec(C)
end
heatmap(string.(collect(1:13)), string.(collect(1:13)),symQ, xlabel=" values of d ", ylabel= " values of k (note: k <= d)")
plot(symQ, xlabel= "n", xticks=1:13, ylabel= "Total Variation")
for i in 1:13 symQ[:,i] /= sum(abs.(symQ[:,i])) end
plot(symQ[:,1], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="1st order")
plot!(symQ[:,2], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="2nd order")
plot!(symQ[:,3], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="3rd order")
plot!(symQ[:,4], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="4th order")
plot!(symQ[:,5], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="5th order")
plot!(symQ[:,6], xlabel= "n", xticks=1:13, ylabel= "Total Variation", label="6th order")


plot()
P[:,1:3]
println(P[:,6])

println(D)
heatmap(string.(judOrigs), string.(judOrigs),D)
#=
274.0   205.0   92.0    127.0   52.0    107.0   39.0
211.0   295.0   126.0   131.0   43.0    159.0   58.0
96.0    106.0   127.0   40.0    20.0    78.0    17.0
119.0   137.0   29.0    109.0   21.0    77.0    22.0
40.0    48.0    14.0    29.0    65.0    51.0    10.0
121.0   179.0   77.0    61.0    42.0    208.0   37.0
35.0    53.0    19.0    17.0    14.0    45.0    93.0
=#
E = D ./ (D+D')
heatmap(string.(judOrigs), string.(judOrigs),E, title="Normalized Density Matrix")

# x is right Perron vector, y is left Perron vector
# note: Perron vector is just the eigenvector where each entry has the same sign
# see Horn P.546 8.6.2
# http://www.cse.zju.edu.cn/eclass/attachments/2015-10/01-1446086008-145421.pdf
x = eigvecs(D)[:,7]
y = eigvecs(D')[:,7]

x*y'

#=
H_0 : π ~ Unif(S_n)
Note: If this is our null hypothesis we could simply rephrase it as: D is symmetric
MLE is
of Prob over S_n is

E(|H_0

=#

D * ones(7)
F = 2* D ./ (D+D')
svd(F)
(D - D') ./2

tr(D*D')

F*F'


for w in waves
    S = Iterators.product([w.panel[c] for c in ]...)
    for τ in permGroup(judOrigs)

waves[3].panel
