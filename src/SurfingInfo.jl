include("SymGrpAndReps.jl")

using JSON: parse

function totVar(A::Array{T,2}) where T <: Number
    if size(A)[1] != size(A)[2] error("A is not square") end
    n = size(A)[1]
    v = 0
    for τ in S(n)
        v += abs(reduce(*,[ A[i,τ[i]] for i in 1:n]) )
    end
    return v
end Defining isordered on orderable collections of numbers
function isordered(X::Array{T}) where T <: Number
    return all([ X[i] < X[i+1] for i in 1:(length(X)-1) ])
end

function isordered(X::Tuple{T}) where T <: Number
    return all([ X[i] < X[i+1] for i in 1:(length(X)-1)])
end

function isordered(X::NTuple{N,T}) where T <: Number where N
    return all([ X[i] < X[i+1] for i in 1:(N-1) ])
end

# Defining isordered on orderable collections of collections
# First given two arrays
function isordered(X::Array{T},Y::Array{T}) where T <: Number
    return all(isordered.( Base.product(X,Y) ) )
end

# Given an array of arrays
function isordered(X::Array{Array{T,1},1}) where T <: Number
    return all( isordered.(Base.product(X...) ))
end

function subtypetree(t, level=1, indent=2)
    level == 1 && println(t)
    for s in subtypes(t)
        println(join(fill(" ", level * indent)) * string(s))
        subtypetree(s, level+1, indent)
    end
end

isoDict = Dict([
    "Australia" => :AUS,
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
    "United States" => :USA
])
ISOs = sort(unique(collect(values(isoDict))))

data = parse(open("../Data/CombinedCountries/CleanAllDataCC.txt", "r"))

waves = []
for wid in keys(data)
    if data[wid]["nJudOrigs"] == 5 & data[wid]["nSubScos"] == 5
        origs = unique(data[wid]["subScoOrig"])
        matchIndicator = (data[wid]["athOrig"] in origs)

        labeledScos = Dict([isoDict[origin] => Float16[] for origin in origs])
        labeledScosBinary = Dict([:Match => Float16[], :NoMatch => Float16[] ])

        origScoPairs = zip(data[wid]["subScoOrig"], data[wid]["subSco"])
        for p in origScoPairs
            # PsuedoCode: push!( array of judge scores from country p[1], score=p[2] )
            push!(labeledScos[ isoDict[p[1]] ], p[2])
            if p[1] == data[wid]["athOrig"]
                push!(labeledScosBinary[:Match], p[2])
            else
                push!(labeledScosBinary[:NoMatch], p[2])
            end
        end

        wave = (
            id=wid,
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
            subScoOrigs= map(x->isoDict[x], data[wid]["subScoOrig"]),
            panelOrigs= Set(map(x->isoDict[x], data[wid]["subScoOrig"])),
            match= matchIndicator
        )

        push!(waves, wave)
    end
end

AthIds = sort(unique(map(x-> x.athId, waves)))
EvtIds = sort(unique(map(x-> x.evtId, waves)))
RndIds = sort(unique(map(x-> x.rndId, waves)))
HeatIds = sort(unique(map(x-> x.heatId, waves)))
WaveIds = sort(map(x-> x.id, waves))