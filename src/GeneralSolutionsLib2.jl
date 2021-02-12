using Plots
gr()
using Distributions
using Statistics
using LinearAlgebra
using JSON
using Iterators

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

isoS = collect(values(isoDict))
e = sort(unique(isoS))
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

w = first(waves)


w.panelBinary

function isordered(X::Array{T},Y::Array{T}) where T <: Number
    if all(map(t-> isless(t[1],t[2]), Base.product(X,Y)))
        return true
    else
        return false
    end
end

function strictOrder(S)
    for i in 1:(length(S)-1)
        if !(S[i] < S[i+1])
            return false
        end
    end
    return true
end

stdBinaryPanels = filter(w->all(length.(values(w.panelBinary)) .> 0), waves)
for country in isoS
    DataPanelBinary = map(w-> w.panelBinary, filter(x-> x.athOrig == country, stdBinaryPanels))
    noOrderBinary = 0
    yesOrderBianry = 0
    DBinary = zeros(2,2)
    for pb in DataPanelBinary
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

A = []
4 = A[:USA, :Bra]

DataPanel = map(x-> x.panel, filter(w->all(length.(values(w.panel)) .> 0), waves))
noOrder = 0
D = zeros(2,2)
for pb in DataPanelBinary
    if isordered(pb[:Match],pb[:NoMatch])
        D += [1 0; 0 1]
    elseif isordered(pb[:NoMatch], pb[:Match])
        D += [0 1; 1 0]
    else
        noOrder += 1
    end
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
    for w in filter(x -> x.panelOrigs == panelComp, waves)
        e = sort(collect(panelComp))
        waveMatrix = zeros(7,7)
        for ord in permGroup(collect(panelComp))
            S = Iterators.product([w.panel[c] for c in ord]...)
            if all(map(strictOrder, S) )
                for (i,ei) in enumerate(e)
                    waveMatrix[ NatToIndex[ei], NatToIndex[ ord[i]] ] = 1
                end
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
T = reduce(+,W)

D = zeros(7,7)
for r in 1:7
    for c in 1:7
        D[r,c] = densityOp[ (judOrigs[r],judOrigs[c]) ]
    end
end
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
