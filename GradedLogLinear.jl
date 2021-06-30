include("DataAndDefns.jl")
L = map(x->x.λ_c, WAVES)

# This is suff. statistic for ~graded~ log linear model!
# A = ⨁̂a = (⨁a¹,…,⨁a^{\#})
A = ⨁(map(l-> hcat( [ sum(e.(pod)) for pod in l ]... ),L))
# 𝔸 = A
𝔸 = hcat([ sum([ A[k][:,j] for k in j:5 ]) for j in 1:5 ]...)

# N = number of judges
N = sum(𝔸)
# Judge Orig margingal is
margJudOrig = 𝔸*ones(5)
# class marginal is
margCls = ones(7)' * 𝔸

# estimated joint distribution assuming independence is
estJoint = margJudOrig * margCls / N
# Chi Squared
sum((𝔸 - estJoint) .^2 ./ estJoint)

ℙJudOrig = margJudOrig/N
ℙCls = margCls/N
ℙ𝔸 = 𝔸/sum(𝔸)


reducedL = map(L) do l
orgs = unique(union(l...))
k = length(orgs)
idx = Dict([c => i for (i,c) in enumerate(sort(orgs)) ])
q = map(cls->sum([e(idx[c],k) for c in cls])/length(cls),l)
reduce(⊗,q)
end;
GL = ⨁(reducedL);
size.(GL)