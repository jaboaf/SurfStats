using Plots; gr()
using Distributions
using Statistics
using LinearAlgebra
using JSON


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

function totVar(A)
    n = size(A)[1]
    if n >= 14 error("Are you sure you want to do that?")
    v = 0
    for τ in permGroup( collect(1:n) )
        v += abs(prod([A[ i , τ[i] ] for i in 1:n]))
    end
    return v
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

data = JSON.parse(open("Data/CombinedCountries/CleanAllDataCC.txt", "r"))
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
                actualSco=data[wid]["actualSco"],
                subScos=data[wid]["subSco"],
                subScoOrigs=map(x->isoDict[x], data[wid]["subScoOrig"]),
                panelOrigs=Set(map(x->isoDict[x], data[wid]["subScoOrig"])),
                match=matchIndicator )

        push!(waves, x)
    end
end


allPanelCompositions = unique( map(x->x.panelOrigs, waves))
judOrigs = sort(collect(∪(allPanelCompositions...)))
permGroup(judOrigs)
NatToIndex = Dict([orig => i for (i,orig) in enumerate(judOrigs) ])
U = 1/7 * ones(7,7)
W = []
totalConsistentPanels = 0
for panelComp in allPanelCompositions
    consistentPanels = 0
    notOnPanel = map(x->NatToIndex[x], setdiff(judOrigs,panelComp) )
    unifM = zeros(7,7)
    unifM[notOnPanel,notOnPanel] .= 1 / length(notOnPanel)
    for w in filter(x -> x.panelOrigs == panelComp, waves)
        e = sort(collect(panelComp))
        waveMatrix = zeros(Int,7,7)
        hasOrder = false
        for ord in permGroup(collect(panelComp))
            S = Iterators.product([w.panel[c] for c in ord]...)
            if all(map(strictOrder, S) )
                for (i,ei) in enumerate(e)
                    waveMatrix[ NatToIndex[ei], NatToIndex[ ord[i]] ] = 1
                end
                waveMatrix += unifM
                consistentPanels +=1
                hasOrder = true
                push!(W, waveMatrix)
            end
        end
        #if hasOrder == false push!(W, U) end
    end
    totalConsistentPanels += consistentPanels
    println(panelComp)
    println(consistentPanels)
end

nWaves = length(W)
origLabels = string.(judOrigs)

totVar(length(W)^-1 * sum(W) - U)

Moments = [ (length(W))^-1 * reduce(+, (W) .^ k ) for k in 1:13 ]
heatmap(origLabels, origLabels, Moments[  5 ] )

CentralMoments = [ (nWaves)^-1 * sum(map(w -> (w - Moments[1]) .^ d, W)) for d in 1:13 ]
heatmap(origLabels, origLabels, CentralMoments[ 5 ] )
