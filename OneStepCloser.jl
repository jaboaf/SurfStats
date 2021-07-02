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
	wv = (
		szn = evtYearToSZN[wave[2]["evtYear"]],
		evt = evtNameToEVT[wave[2]["evtYear"]],
		rnd = wave[2]["rnd"],
		ht = wave[2]["heat"],
		athName = wave[2]["athName"],
	subScoOrig = map(x-> isoDict[x], wave[2]["subScoOrig"])
	subSco = map(x->round(x,digits=1),wave[2]["subSco"])
	athorig = isoDict[ wave[2]["athOrig"] ]

	origs = unique(subScoOrig)
	scos = Tuple(sort(unique(subSco)))

	part_c = map(s-> subScoOrig[ findall(==(s),subSco) ], scos)
	part_b = map(X-> X .== athorig ,part_c)

	mult_c = map(c->count(==(c),subScoOrig), C)
	mult_b = map(b->count(==(b),subScoOrig .== athorig ), (0,1) )

