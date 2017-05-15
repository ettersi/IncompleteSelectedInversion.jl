module IncompleteSelectedInversion

#immutable UnsafeArrayWrapper{T,N} <: AbstractArray{T,N}
#    data::Ptr{T}
#    dims::NTuple{N,Int}
#    IsbitsArrayWrapper(a::AbstractArray{T,N})
#end


#=
 Sorted subset of {1,...,n}
 next[n+1] = first entry in the set
 If i is an element of the set, then next[i] is the 
 smallest j > i in the set, or n+1 if no such j exists.
=#
immutable SortedIntSet
    next::Vector{Int} 
    SortedIntSet(n) = new(Vector{Int}(n+1))
end

#= 
 Iteration
=#
Base.start(s::SortedIntSet) = length(s.next)
Base.next(s::SortedIntSet,p) = s.next[p],s.next[p]
Base.done(s::SortedIntSet,p) = s.next[p] == length(s.next)

function init!(s::SortedIntSet,i)
    next = s.next
    n = length(next)-1
    next[n+1] = i
    next[i] = n+1
    return s
end

#=
 Insert i into the set. p points to an element before i in the linked list. 
 Return whether the element was newly inserted. 
=#
Base.insert!(s::SortedIntSet,i) = insert!(s,i,length(s.next))
function Base.insert!(s::SortedIntSet,i,p)
    next = s.next
    n = length(next)-1

    if !(1 <= i <= n); throw(BoundsError()); end
    if !(p == n+1 || 1 <= p <= i); throw(BoundsError()); end
    @inbounds begin
        while next[p] < i
            p = next[p]
        end
        if next[p] == i
            return false
        end
        next[p],next[i] = i,next[p]
        return true
    end
end

abstract AbstractSparsePI{Ti}
ncols(A::AbstractSparsePI) = length(A.p)-1
idxtype{Ti}(A::AbstractSparsePI{Ti}) = Ti


immutable SparsePI{Ti} <: AbstractSparsePI{Ti}
    p::Vector{Ti}
    i::Vector{Ti}
end

immutable SparsePIL{Ti} <: AbstractSparsePI{Ti}
    p::Vector{Ti}
    i::Vector{Ti}
    l::Vector{Int}
end
function (::Type{SparsePIL{Ti}}){Ti}(n)
    p = Vector{Ti}(n+1)
    p[1] = 1
    i = Vector{Ti}(0)
    l = Vector{Ti}(0)
    return SparsePIL{Ti}(p,i,l)
end

immutable SparsePIL_Column
    i::SortedIntSet
    l::Vector{Int}
    c::Int
end
function (::Type{SparsePIL_Column})(n,c)
    i = SortedIntSet(n)
    l = Vector{Int}(n)
    return SparsePIL_Column(i,l,c)
end

immutable SparsePIV{Ti,Tv}
    p::Vector{Ti}
    i::Vector{Ti}
    v::Vector{Tv}
end

function init!(Fj::SparsePIL_Column,j,A)
    init!(Fj.i,j)
    Fj.l[j] = 0 
    lasti = j
    for p in A.p[j]:A.p[j+1]-1
        i = A.i[p]
        if i <= j; continue; end
        insert!(Fj.i,i,lasti)
        Fj.l[i] = 0 
        lasti = i
    end
    return nothing
end
function update!(Fj::SparsePIL_Column,F::SparsePIL,pvals)
    n = ncols(F)
    lkj = F.l[first(pvals)]
    if lkj >= Fj.c; return; end
    lasti = n+1
    for p in pvals
        i = F.i[p]
        lik = F.l[p]
        Flij = lik + lkj + 1
        if Flij <= Fj.c
            if insert!(Fj.i,i,lasti)
                Fj.l[i] = Flij
            else
                Fj.l[i] = min(Fj.l[i],Flij)
            end
            lasti = i
        end
    end
end
function Base.push!(F::SparsePIL,j,Fj::SparsePIL_Column)
    for i in Fj.i
        push!(F.i,i)
        push!(F.l,Fj.l[i])
    end
    F.p[j+1] = length(F.i)+1
    return nothing
end

immutable Iteration_jkp{Structure_A}
    A::Structure_A
    nextp::Vector{Int}
    nextk::Vector{Int}
end
immutable Iteration_kp{Ti}
    jkp::Iteration_jkp{Ti}
    j::Int
end

function iterate_jkp(A::AbstractSparsePI)
    n = ncols(A)
    nextp = Vector{Int}(n)
    nextk = Vector{Int}(n)
    fill!(nextk,n+1)
    fill!(nextp,0)
    Iteration_jkp(A,nextp,nextk)
end
Base.start(jkp::Iteration_jkp) = 0
function Base.next(jkp::Iteration_jkp,j) 
    A = jkp.A
    nextp = jkp.nextp
    nextk = jkp.nextk

    if j > 0
        for p in A.p[j]:A.p[j+1]-1
            i = A.i[p]
            if i > j
                nextp[j] = p
                nextk[i],nextk[j] = j,nextk[i]
                break
            end
        end
    end
    j += 1
    return (j,Iteration_kp(jkp,j)),j
end
Base.done(jkp::Iteration_jkp,j) = j == ncols(jkp.A)

Base.start(kp::Iteration_kp) = kp.jkp.nextk[kp.j]
function Base.next(kp::Iteration_kp,k)
    jkp = kp.jkp
    A = jkp.A
    nextp = jkp.nextp
    nextk = jkp.nextk

    pp = nextp[k]
    kk = nextk[k]
    nextp[k] += 1
    if nextp[k] < A.p[k+1]
        i = A.i[nextp[k]]
        nextk[i],nextk[k] = k,nextk[i]
    end
    return (k,pp:A.p[k+1]-1),kk
end
Base.done(kp::Iteration_kp,k) = k > kp.j


function symbolic(A,c)
    Ti = idxtype(A)
    n = ncols(A)

    F = SparsePIL{Ti}(n)
    Fj = SparsePIL_Column(n,c)
    for (j,kvals) in iterate_jkp(F)
        init!(Fj,j,A)
        for (k,pvals) in kvals
            update!(Fj,F,pvals)
        end
        push!(F,j,Fj)
    end
    return F
end

end # module
