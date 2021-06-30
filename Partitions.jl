

PartData = map(x->sort(collect(length.(x.λ_c)),rev=true),Tups);
Parts = sort(unique(PartData),rev=true)
PartCounts = [count(==(p),PartData) for p in Parts]
barplot(Parts,PartCounts,title="Empirical Distribution of Unordered Panel Partitions")

OrdPartData = map(x->length.(x.λ_c),Tups);
OrdParts = sort(unique(OrdPartData),rev=true)
OrdPartCounts = [count(==(p),OrdPartData) for p in OrdParts]
barplot(OrdParts,OrdPartCounts,title="Empirical Distribution of Ordered Panel Partitions")

# [2,1,1,1]
barplot(OrdParts[[2,3,5,9]],OrdPartCounts[[2,3,5,9]],title="Distribution of Orders of [1,1,3]")
#[2,2,1]
barplot(OrdParts[[6,10,11]],OrdPartCounts[[6,10,11]])
#[3,1,1]
barplot(OrdParts[[4,7,13]],OrdPartCounts[[4,7,13]])
#[3,2]
barplot(OrdParts[[12,14]],OrdPartCounts[[12,14]])
#[4,1]
barplot(OrdParts[[8,15]],OrdPartCounts[[8,15]])

for r in R
	map()

