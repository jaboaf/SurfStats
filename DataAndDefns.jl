include("SymGrpAndReps.jl")
using JSON

# 3 instances of SZN
@enum SZN WCT17=1 WCT18 WCT19
evtYearToSZN = Dict([ "2017"=>WCT17,"2018"=>WCT18,"2019"=>WCT19 ])
# 12 instances of EVT
@enum EVT BaliPro=1 BellsBeach Fiji FrancePro GoldCoast JBayOpen MargaretRiver PenichePro PipeMasters RioPro Teahupoo Trestles
evtNameToEVT = Dict([
	"Bali Pro" => BaliPro, 
	"Bells Beach" => BellsBeach,
	"Fiji" => Fiji,
	"France" => FrancePro,
	"Gold Coast" => GoldCoast,
	"J-Bay Open" => JBayOpen,
	"Margaret River" => MargaretRiver,
	"Peniche Pro" => PenichePro,
	"Pipe Masters" => PipeMasters,
	"Rio Pro" => RioPro,
	"Tahiti" => Teahupoo,
	"Trestles" => Trestles
])

const RND = [i for i in 1:7]
const CC = false
if CC
	data = JSON.parse( open("data/CleanAllDataCC.txt", "r") )
	@enum ORIG AUS=1 BRA ESP FRA PRT USA ZAF FJI IDN ITA JPN NZL
	const C = collect(instances(ORIG)[1:7])
	#= CTRY Defined by ISO 3166-1
	Using 2 digit codes
	"-1"
	=#
	isoDict = Dict([
		"Australia" => AUS,
		"Brazil" => BRA,
		"Basque Country" => ESP,
		"Spain" => ESP,
		"France" => FRA,
		"Portugal" => PRT,
		"United States" => USA,
		"South Africa" => ZAF,
		"Fiji" => FJI,
		"Indonesia" => IDN,
		"Italy" => ITA,
		"Japan" => JPN,
		"New Zealand" => NZL
	])
else
	data = JSON.parse( open("data/CleanAllDataNC.txt", "r") )
	@enum ORIG AUS=1 BAS BRA FRA HAW PRT PYT USA ZAF ESP FJI IDN ITA JPN NZL
	const C = instances(ORIG)[1:9]
	#= ORIG Defined by what WSL had. 
	Using 3 digit codes from ISO 3166-1, plus BAS for Basque Country.
	"-1"
	=#
	isoDict = Dict([
		"Australia" => AUS,
		"Basque Country" => BAS,
		"Brazil" => BRA,
		"France" => FRA,
		"French Polynesia" => PYT,
		"Hawaii" => HAW,
		"Portugal" => PRT,
		"United States" => USA,
		"South Africa" => ZAF,
		"Spain" => ESP,
		"Fiji" => FJI,
		"Indonesia"=>IDN,
		"Italy"=>ITA,
		"Japan"=>JPN,
		"New Zealand"=>NZL
	])
end

filter!(wave-> wave[2]["subScoOrigDefect"]==false,data);

WIDs= sort(collect(keys(data)));
EvtIds = unique(map(x->data[x]["evtId"],WIDs));
HeatIds = unique(map(x->data[x]["heatId"],WIDs));
maxRnd = Dict()
for id in EvtIds
	thisEvt = filter(w->data[w]["evtId"]==id, WIDs)
	maxRnd[id] = maximum(map(w->Base.parse(Int,data[w]["rnd"]), thisEvt))
end

WAVES = []
for wave in data
	subScoOrig = map(x-> isoDict[x], wave[2]["subScoOrig"])
	subSco = map(x->round(x,digits=1),wave[2]["subSco"])
	athorig = isoDict[ wave[2]["athOrig"] ]

	origs = unique(subScoOrig)
	scos = Tuple(sort(unique(subSco)))

	part_c = map(s-> subScoOrig[ findall(==(s),subSco) ], scos)
	part_b = map(X-> X .== athorig ,part_c)

	mult_c = map(c->count(==(c),subScoOrig), C)
	mult_b = map(b->count(==(b),subScoOrig .== athorig ), (0,1) )

	wv = (
		id = wave[1],
		szn = evtYearToSZN[ wave[2]["evtYear"] ],
		evtOrig = isoDict[ wave[2]["evtOrig"] ],
		evt = evtNameToEVT[ wave[2]["evtName"] ],
		rnd = maxRnd[wave[2]["evtId"]]-Base.parse(Int,wave[2]["rnd"])+1,
		ht = Base.parse(Int, wave[2]["heat"]),
		athOrig = isoDict[ wave[2]["athOrig"] ],
		I_match = athorig in subScoOrig,
		m_c = mult_c,
		m_b = mult_b,
		λ_c = part_c,
		λ_b = part_b
	)
	push!(WAVES,wv)
end
sort(WAVES);


⊗(A::Array{T},B::Array{T}) where T<: Number = prod.(Base.product(A,B))
⊗(V::Vararg{Array{T}}) where T<: Number =reduce(⊗,V) 

⊗(a::NTuple{T},b::NTuple{T}) where T  = (a...,b...)
⊗(a::NTuple{N,T},b::NTuple{M,T}) where {T,N,M}  = (a...,b...)

⊕(V::Vararg{Array{T,N} where N}) where T = [ sum(filter(x->size(x)==d,V)) for d in unique(size.(V)) ] 
⊕(V::Vararg{Array{T,N}}) where {T,N} = [ sum(filter(x->size(x)==d,V)) for d in unique(size.(V)) ] 


# maybe
⨁(A::Array{Array{T,N} where N,1}) where T <: Number = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨁(A::Array{Array{T},1}) where T<:Number = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨁(A::Array{Array{T,N},1}) where {T<:Number,N} = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨂(A::Array{Array{T,N} where N,1}) where T <: Number = reduce(⊗,A)


E(X::Array,i::K) where {K<:Integer} =dropdims( sum(X,dims=setdiff(1:ndims(X),i)),dims=tuple(setdiff(1:ndims(X),i)...) )
E(X::Array,I::NTuple) =dropdims(sum(X,dims=setdiff(1:ndims(X),I)),dims=tuple(setdiff(1:ndims(X),I)...))

μₐ(D::Array{Array{T,N} where N}) where T<:Number = map(x->sum(abs.(S(x))),⨁(D))
μₛ(D::Array{Array{T,N} where N}) where T<:Number = map(x->sum(abs.(S(x))),⨁(D))

#cov(X::Array,i::K,j::K) where {K<:Integer} = E(X,(i,j))-E(X,i)⊗E(X,j)
#cov(X::Array,I::NTuple,j::K) where {K<:Integer} = E(X,(I...,j))-E(X,I)⊗E(X,j)
#cov(X::Array,i::K,J::NTuple) where {K<:Integer} = E(X,(i,J...))-E(X,i)⊗E(X,J)
#cov(X::Array,I::NTuple,J::NTuple) =E(X,(I...,J...))-E(I)⊗E(J)

e_orig(c::ORIG) = vcat(zeros(Int(c)-1),1,zeros(length(instances(ORIG))-Int(c)))
e(c::ORIG) = vcat(zeros(Int(c)-1),1,zeros(length(C)-Int(c)))
e(b::Bool) = [ Int(b==0); Int(b==1)] 
e(i::Integer,n::Integer) = vcat(zeros(i-1),1,zeros(n-i)) 
function embedd(t::Tuple{T,Vararg{T,N} where N} where T <:Union{Array{ORIG,1},BitArray{1}};strong=false, countinblock=false, minimal=false, normalize=false)
	if minimal==true # then embedd in VS with as many dims are needed
		orgs = unique(union(t...))
		k = length(orgs)
		idx = Dict([c => i for (i,c) in enumerate(sort(orgs)) ])
		if countinblock==true
			q = map(cls->[ count(==(c),cls)*e(idx[c],k) for c in unique(cls)],t)
		else
			q = map(cls->map(c->e(idx[c],k),cls),t)
		end
	else # embedd in VS with e_Int(label)
		if countinblock==true
			q = map(cls->[count(==(c),cls)*e(c) for c in unique(cls)],t)
		else
			q = map(cls->e.(cls),t)
		end
	end
	if normalize==true
		q = map(pod -> pod ./ sum(sum(pod)), q)
	end
	if strong==true # then do sum within pod
		q = map(pod->sum(pod),q)
		return ⊗(q...)
	else
		blocks = map(pod->S(⊗(pod...)),q)
		return ⊗(blocks...)
	end
end















