function Sym(A::Union{Set,Array})
    if length(A) >= 14
        error("Slow down... be careful, n >= 14, and |S_n| = $(factorial(big(length(A)))) ")
    end
    unique!(A)
    Perms = Array{eltype(A),1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(Perms, head)
        end
    end
    continuePerm(eltype(A)[], A)
    return Perms
end

function Sym(n::Integer)
    if n >= 14
        error("Slow down... be careful, n >= 14, and |S_n| = $(factorial(big(n))) ")
    end
    A = collect(Int8, 1:n)
    Perms = Array{Int8,1}[]
    function continuePerm(head,tail)
        if length(tail) > 0
            for t in tail
                newHead = union(head, [t])
                newTail = setdiff(tail, [t])
                continuePerm(newHead, newTail)
            end
        else
            push!(Perms, head)
        end
    end
    continuePerm(Int8[], A)
    return Perms
end

function Rep(t::Array)
    M = zeros(Bool, length(t),length(t))
    for i in 1:length(t)
        M[ t[i] , i ] = 1
    end
    return M
end


function Rep(t::Array, Res::Array)
    if !(issubset(Res, t))
        error("Restriction is not contained in the domain of the perm, t.")
    end
    M = zeros(Bool, length(t), length(t))
    for j in Res
        M[ t[j] , j ] = 1
    end
    return M
end

function Rep(t::Array, Ind::Array)
    if !(issubset(t, Ind))
        error("The image of the perm t is not contained in the Induced space")
    end
    M = zeros(Bool, length(Ind), length(Ind))
    for j in sort(t)
        M[ t[j] , j ] = 1
    end
    return M
end