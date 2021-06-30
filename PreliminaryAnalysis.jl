include("DataAndDefns.jl")
using Statistics
using GRUtils

#' The straight forward approach is to test the differences in means between judges with the same nationality as the surfer and those with a different nationality. I.e. Test
#' H₀: Mean(Match Scores) - Mean(No Match Scores) = 0
#' H₁: Mean(Match Scores) - Mean(No Match Scores) != 0
ttest(X) = mean(X) / (std(X)/sqrt(length(X)))

matchWaves = filter(x->x.I_match==true, WAVES)
Diffs = Dict([c => Float64[] for c in C])
for m in matchWaves
	subscoorig = map(x-> isoDict[x], data[m.id]["subScoOrig"])
	subsco = map(x->round(x,digits=1),data[m.id]["subSco"])
	I = findall(==(m.athOrig),subscoorig)
	d = mean( subsco[I] ) - mean( subsco[setdiff(1:5,I)] )
	push!(Diffs[m.athOrig], d)
end

# We may carry out this process for each country:
function diffInMeansDistribution(c::ORIG)
	N = length(Diffs[c])
	t = round(ttest(Diffs[c]), digits=6)
	m = round(mean(Diffs[c]), digits=6)
	savefig(
		"visuals/ttests/CC$(CC)/DistOfDiffInMeansFor$(c).png",
		histogram(Diffs[c],title="Difference in Means Given Origin = $c (N=$N, μ̂=$m, t=$t)")
	)
	return t
end
for c in C if length(Diffs[c])>0 diffInMeansDistribution(c) end end

Figure()
Ctry = []
N = []
M = []
T = []
for c in C
	if length(Diffs[c]) > 0
		push!(Ctry, " $c ")
		n = length(Diffs[c])
		m = round(mean(Diffs[c]),digits=3)
		t = round(ttest(Diffs[c]), digits=3)
		push!(N," $n ")
		push!(M," $m ")
		push!(T," $t ")
		plot(sort(Diffs[c]),LinRange(0,1,length(Diffs[c])), label="$c ",hold=true,xlim=(-1.5,1.5))
	end
end
if CC title("Empirical CDF of Difference in Means by Country")
else title("Empirical CDF of Difference in Means by Origin")
end
println(join(Ctry, "&"), " \\\\ \\hline")
println(join(N, "&"), " \\\\ \\hline")
println(join(M, "&"), " \\\\ \\hline")
println(join(T, "&"), " \\\\ \\hline")
legend(location=4)
savefig("visuals/ttests/CC$(CC)/EmpDiffInMeansCDF.png",gcf())
