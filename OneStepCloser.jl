include("DataAndDefns.jl")

EXACT = []
for id in WIDs
	subSco = map(x->round(x,digits=1),data[id]["subSco"])
	subScoOrig = map(x-> isoDict[x], data[id]["subScoOrig"])
	scos = sort(unique(subSco)) 
	part = map(s-> subScoOrig[ findall(==(s),subSco) ], scos)
	for (i,s) in enumerate(scos)
		subsc = (
			szn = evtYearToSZN[data[id]["evtYear"]],
			evt = evtNameToEVT[data[id]["evtName"]],
			rnd = data[id]["rnd"],
			ht = data[id]["heat"],
			athorig = isoDict[data[id]["athOrig"]],
			athname = data[id]["athName"],
			sco = s,
			jud = part[i]
		)
		push!( EXACT, subsc)
	end
end
L = map(x->x.λ_c,WAVES);
K = vcat(collect.(L)...);

# szn in (WCT17, WCT18, WCT19)
# evt in instances(EVT)
# rnd in RND = (1,\dots,8)
# ht in HT = (1,\dots,16)
# athorig in instances(ORIG)
# athname in ATHNAME 
ATHNAME = Tuple(sort(unique(map(x->x.athName,WAVES))))
# s in SCORE = (1,\dots,100)
SCORE = Tuple(1:100)

D = Array{Array{Rational,N} where N}(undef,(3,12,8,16,15,88,100))


𝕐 = e_G.(L);
𝐒𝕐 = 𝐒.(𝕐);
𝕖 = 𝕐 - 𝐒𝕐