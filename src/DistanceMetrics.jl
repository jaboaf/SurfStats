import JSON

data = Dict()
open("Data/CombinedCountries/CleanAllDataCC.txt", "r") do f
    global data
    data = JSON.parse(f)  # parse and transform data
end

# some notes
data :: Dict{String, Any}
data[somekey] :: Dict{String, Any}


origs = map( x -> x["subScoOrig"], values(data) )
filter!( x -> !("-1" in x), origs)
unique(origs) # 574

unique(map( x -> sort(collect(values(countmap(x)))), origs ))
# Sorted counts of countries

countmap(map( x -> sort(collect(values(countmap(x)))), origs ))
#=
[1, 1, 1, 1, 1] => 3234
[1, 1, 1, 2]    => 8284
[1, 2, 2]       => 3013
[1, 1, 3]       => 1216
[2, 3]          => 225
[1, 4]          => 18
=#

countmap(map( x -> collect(values(countmap(x))), origs ))
#=
DISCLAIMER: think about what we are actually looking at here, (values of countmap)
What do we actually want?
   first country seen -> a
   second country seen -> b
   ...
   5th country seen -> e
   (if a country was already seen, we don't "allocate" a new letter )
   then replace each country in list by the corresponding letter
   return all unique (label-blind) lists

   ... Maybe we should be testing inverse permutation of rankings

[1, 1, 1, 1, 1]   => 3234

[2, 1, 1, 1]      => 835
[1, 2, 1, 1]      => 3679
[1, 1, 2, 1]      => 2857
[1, 1, 1, 2]      => 913

[1, 2, 2]         => 1180
[2, 1, 2]         => 217
[2, 2, 1]         => 1616

[3, 1, 1]         => 608
[1, 3, 1]         => 274
[1, 1, 3]         => 334

[4, 1]            => 18
[3, 2]            => 225

=#

filter( x -> !("-1" in x[2]["subScoOrig"] )

end, get.( values(data), "subScoOrig" ) )

# Kendall’s Tau
function KendallTau(τ₁::T, τ₂::T where T::Array)::Float64
   for ls


function KendallTau(τ₁::T, τ₂::T where T::Array{Set{Any}})::Float64
   for

distance between two lists {\displaystyle \tau _{1}}\tau _{1} and {\displaystyle \tau _{2}}\tau _{2} is

   {\displaystyle K(_{1},\tau _{2})=|\{(i,j):i<j,(\tau _{1}(i)<\tau _{1}(j)\wedge \tau _{2}(i)>\tau _{2}(j))\vee (\tau _{1}(i)>\tau _{1}(j)\wedge \tau _{2}(i)<\tau _{2}(j))\}|.}{\displaystyle K(\tau _{1},\tau _{2})=|\{(i,j):i<j,(\tau _{1}(i)<\tau _{1}(j)\wedge \tau _{2}(i)>\tau _{2}(j))\vee (\tau _{1}(i)>\tau _{1}(j)\wedge \tau _{2}(i)<\tau _{2}(j))\}|.}
   where
