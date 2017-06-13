module IncompleteSelectedInversion

#=
 Sorted subset of {1,...,n}.
 Used to represent the row indices of a single column in F.
=#
immutable SortedIntSet
    next::Vector{Int} 
    #=
     next[n+1] = first entry in the set
     If i is an element of the set, then next[i] is the 
     smallest j > i in the set, or n+1 if no such j exists.
    =#
    SortedIntSet(n) = new(Vector{Int}(n+1))
end

Base.start(s::SortedIntSet) = length(s.next)
Base.@propagate_inbounds Base.next(s::SortedIntSet,p) = s.next[p],s.next[p]
Base.@propagate_inbounds Base.done(s::SortedIntSet,p) = s.next[p] == length(s.next)

Base.@propagate_inbounds function init!(s::SortedIntSet,i)
    next = s.next
    n = length(next)-1

    @boundscheck begin
        @assert 1 <= i <= n
    end

    @inbounds begin
        next[n+1] = i
        next[i] = n+1
        return s
    end
end

Base.@propagate_inbounds function Base.insert!(s::SortedIntSet,i,p)
    next = s.next
    n = length(next)-1

    @boundscheck begin
        @assert 1 <= i <= n
        @assert p == n+1 || 1 <= p <= i
    end
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



#=
 Iterate through a sparse matrix in the order required
 by the LDL^T factorisation. 
=#

iterate_jkp(Ap,Ai) = Iteration_jkp(Ap,Ai)

immutable Iteration_jkp{Ti}
    Ap::Vector{Ti}
    Ai::Vector{Ti}
    nextk::Vector{Int}
    nextp::Vector{Int}
end
immutable Iteration_kp{Ti}
    Ap::Vector{Ti}
    Ai::Vector{Ti}
    nextk::Vector{Int}
    nextp::Vector{Int}
    j::Int
end

function Iteration_jkp(Ap,Ai)
    n = length(Ap)-1
    nextk = Vector{Int}(n)
    fill!(nextk,n+1)
    nextp = Vector{Int}(n)
    return Iteration_jkp(Ap,Ai,nextk,nextp)
end
Base.start(jkp::Iteration_jkp) = 0
Base.@propagate_inbounds function Base.next(jkp::Iteration_jkp,j) 
    Ap = jkp.Ap
    Ai = jkp.Ai
    nextp = jkp.nextp
    nextk = jkp.nextk

    if j > 0
        for p in Ap[j]:Ap[j+1]-1
            i = Ai[p]
            if i > j
                nextp[j] = p
                nextk[i],nextk[j] = j,nextk[i]
                break
            end
        end
    end
    j += 1
    return (j,Iteration_kp(Ap,Ai,nextk,nextp,j)),j
end
Base.done(jkp::Iteration_jkp,j) = j == length(jkp.Ap)-1

Base.@propagate_inbounds Base.start(kp::Iteration_kp) = kp.nextk[kp.j]
Base.@propagate_inbounds function Base.next(kp::Iteration_kp,k)
    Ap = kp.Ap
    Ai = kp.Ai
    nextp = kp.nextp
    nextk = kp.nextk

    pp = nextp[k]
    kk = nextk[k]
    nextp[k] += 1
    if nextp[k] < Ap[k+1]
        i = Ai[nextp[k]]
        nextk[i],nextk[k] = k,nextk[i]
    end
    return (k,pp:Ap[k+1]-1),kk
end
Base.done(kp::Iteration_kp,k) = k > kp.j


function checkmat(Ap,Ai)
    n = length(Ap)-1
    for j = 1:n
        @assert 1 <= Ai[Ap[j]] <= n
        for p in Ap[j]+1:Ap[j+1]-1
            @assert Ai[p-1] < Ai[p] <= n
        end
    end
end
function checkmat(Ap,Ai,Ax,Ay...) 
    @assert length(Ax) >= Ap[end]-1
    checkmat(Ap,Ai,Ay...)
end


function symbolic(Ap,Ai,c)
    checkmat(Ap,Ai)

    @inbounds begin
        Ti = eltype(Ap)
        n = length(Ap)-1

        # Return variables
        Fp = Vector{Ti}(n+1); Fp[1] = 1
        Fi = Vector{Ti}(0)
        Fl = Vector{Ti}(0)

        # Workspace for a single column
        Fji = SortedIntSet(n)
        next = Fji.next
        Fjl = Vector{Int}(n)

        # Main algorithm
        for (j,kvals) in iterate_jkp(Fp,Fi)
            # Initialise column
            init!(Fji,j)
            Fjl[j] = 0 
            lasti = j
            for p in Ap[j]:Ap[j+1]-1
                i = Ai[p]
                if i <= j; continue; end
                insert!(Fji,i,lasti)
                Fjl[i] = 0 
                lasti = i
            end

            # Pull updates into L[j:n,j]
            for (k,pvals) in kvals
                lkj = Fl[first(pvals)]
                if lkj >= c; continue; end
                lasti = n+1
                for p in pvals
                    i = Fi[p]
                    lik = Fl[p]
                    Flij = lik + lkj + 1
                    if Flij <= c
                        if insert!(Fji,i,lasti)
                            Fjl[i] = Flij
                        else
                            Fjl[i] = min(Fjl[i],Flij)
                        end
                        lasti = i
                    end
                end
            end

            # Copy temporary column into F
            for i in Fji
                push!(Fi,i)
                push!(Fl,Fjl[i])
            end
            Fp[j+1] = length(Fi)+1
        end
        return Fp,Fi,Fl
    end
end


function numeric(Ap,Ai,Av,Fp,Fi)
    checkmat(Ap,Ai,Av)
    checkmat(Fp,Fi)

    @inbounds begin
        Ti = eltype(Ap)
        Tv = eltype(Av)
        n = length(Ap)-1

        # Return variables
        Fv = Vector{Tv}(length(Fi))

        # Workspace for a single column
        Fjv = Vector{Tv}(n)

        # Main algorithm
        for (j,kvals) in iterate_jkp(Fp,Fi)
            # Initialise column
            for p in Fp[j]:Fp[j+1]-1
                Fjv[Fi[p]] = zero(Tv)
            end
            for p in Ap[j]:Ap[j+1]-1
                Fjv[Ai[p]] = Av[p]
            end

            # Pull updates into L[j:n,j]
            for (k,pvals) in kvals
                f = Fv[Fp[k]]*Fv[first(pvals)]
                for p in pvals
                    # We compute a few dropped fill-ins here. It turns out computing 
                    # and discarding is faster than introducing a branch. 
                    Fjv[Fi[p]] -= Fv[p]*f
                end
            end

            # Copy temporary column into F
            d = Fjv[j]
            Fv[Fp[j]] = d
            for p in Fp[j]+1:Fp[j+1]-1
                Fv[p] = Fjv[Fi[p]]/d
            end
        end
        return Fv
    end
end


function selinv_jki(Fp,Fi,Fv)
    checkmat(Fp,Fi,Fv)

    @inbounds begin
        Ti = eltype(Fp)
        Tv = eltype(Fv)
        n = length(Fp)-1

        # Return variables
        Bv = Vector{Tv}(length(Fi))

        # Workspace for a single column
        Fjv = Vector{Tv}(n)
        Bjv = Vector{Tv}(n)

        # Main algorithm
        for j in reverse(1:n)
            # Initialise column
            for p in Fp[j]+1:Fp[j+1]-1
                Fjv[Fi[p]] = Fv[p]
                Bjv[Fi[p]] = zero(Tv)
            end

            # Pull updates into B[j+1:n,j]
            for p in Fp[j]+1:Fp[j+1]-1
                k = Fi[p]
                Fkj = Fjv[k]
                Bjv[k] -= Bv[Fp[k]]*Fkj
                for p in Fp[k]+1:Fp[k+1]-1
                    i = Fi[p]
                    Fij = Fjv[i]
                    Bjv[i] -= Bv[p] *Fkj
                    Bjv[k] -= Bv[p]'*Fij
                end
            end

            # Copy temporary column into B
            for p in Fp[j]+1:Fp[j+1]-1
                Fjv[Fi[p]] = zero(Tv)
                Bv[p] = Bjv[Fi[p]]
            end

            # Deal with diagonal
            d = inv(Fv[Fp[j]])
            for p in Fp[j]+1:Fp[j+1]-1
                d -= Bv[p]*Fv[p]
            end
            Bv[Fp[j]] = d
        end
        return Bv
    end
end



#=
 Utility functions
=#

function dropfillin(Fp,Fi,Fl,c)
    Ti = eltype(Fp)
    n = length(Fp)-1
    F̃p = Vector{Ti}(n+1)
    F̃p[1] = 1
    F̃i = Vector{Ti}(0)
    F̃l = Vector{Ti}(0)
    for j = 1:n
        for p in Fp[j]:Fp[j+1]-1
            if Fl[p] <= c
                push!(F̃i,Fi[p])
                push!(F̃l,Fl[p])
            end
        end
        F̃p[j+1] = length(F̃i)+1
    end
    return F̃p,F̃i,F̃l
end

unpacksparse(A) = A.colptr,A.rowval,A.nzval
function packsparse(p,i,v)
    n = length(p)-1
    return SparseMatrixCSC(n,n,p,i,v)
end

end # module
