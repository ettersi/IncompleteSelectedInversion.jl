module IncompleteSelectedInversion

# Sorted subset of {1,...,n}
# next[n+1] = first entry in the set
# If i is an element of the set, then next[i] is the 
# smallest j > i in the set, or n+1 if no such j exists.
immutable SortedIntSet
    next::Vector{Int} 
    SortedIntSet(n) = new(Vector{Int}(n+1))
end

# Iteration
Base.start(s::SortedIntSet) = length(s.next)
Base.next(s::SortedIntSet,p) = s.next[p],s.next[p]
Base.done(s::SortedIntSet,p) = s.next[p] == length(s.next)
function empty!(s::SortedIntSet) 
    next = s.next
    n = length(next)-1
    next[n+1] = n+1
    return s
end

# Insert i into the set. p points to an element before i in the linked list. 
# Return whether the element was newly inserted. 
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

function symbolic(Ap,Ai,c)
    n = length(Ap)-1
    Ti = eltype(Ap)


    # ----------------
    # Output variables
    # ----------------

    Fp = Vector{Ti}(n+1) 
    Fi = Vector{Ti}(0)
    # Fi[Fp[j]:Fp[j+1]-1] = find(F[j:end,j])

    Fl = Vector{Int}(0) 
    # Fl[p] = level of fill of F[i,j] where i = Fi[p] and p in Fp[j]:Fp[j+1] 

    Fq = Vector{Ti}(n+1) 
    Fk = Vector{Ti}(0)
    # Fk[Fq[i]:Fq[i+1]-1] = find(F[i,1:i-1]), in arbitrary order


    # ---------
    # Workspace
    # ---------

    Fij = SortedIntSet(n) 
    # Fij = find(F[j:n,j])

    Flj = Vector{Int}(n) 
    # Flj[i] = level-of-fill of F[i,j], or garbage if !(i in Fij)

    nextp = Vector{Int}(n)
    nextk = Vector{Int}(n)
    fill!(nextk,n+1)

    Fp[1] = 1
    Fq[1] = 1
    for j = 1:n
        # Initialise Fij = find(A[j:n,j]), Flj[Fij] = 0
        empty!(Fij)
        insert!(Fij,j,n+1); lasti = j
        Flj[j] = 0 
        for i in @view Ai[Ap[j]:Ap[j+1]-1]
            if i <= j; continue; end
            insert!(Fij,i,lasti); lasti = i
            Flj[i] = 0 
        end

        k = nextk[j]
        while k < j
            if Fl[nextp[k]] < c
                lasti = n+1
                for p in nextp[k]:Fp[k+1]-1
                    i = Fi[p]
                    l = Fl[p]
                    Flij = l + Fl[nextp[k]] + 1
                    if Flij <= c
                        # TODO: try saving time by eliminating branch
                        if insert!(Fij,i,lasti)
                            Flj[i] = Flij
                        else
                            Flj[i] = min(Flj[i],Flij)
                        end
                        lasti = i
                    end
                end
            end

            push!(Fk,k)
            newk = nextk[k]
            nextp[k] += 1
            if nextp[k] < Fp[k+1]
                i = Fi[nextp[k]]
                nextk[i],nextk[k] = k,nextk[i]
            end
            k = newk
        end

        # Store Fij and Flj into permanent storage
        for i in Fij
            push!(Fi,i)
            push!(Fl,Flj[i])
        end
        Fp[j+1] = length(Fi)+1
        Fq[j+1] = length(Fk)+1

        for p in Fp[j]:Fp[j+1]-1
            i = Fi[p]
            if i > j
                nextp[j] = p
                nextk[i],nextk[j] = j,nextk[i]
                break
            end
        end
    end
    return Fp,Fi,Fl,Fq,Fk
end

function numeric(Ap,Ai,Av,Fp,Fi,Fq,Fk)
    n = length(Ap)-1
    Ti = eltype(Ap)
    Tv = eltype(Av)

    # ----------------
    # Output variables
    # ----------------

    Fv = Vector{Tv}(length(Fi))

    # ---------
    # Workspace
    # ---------

    Fvj = Vector{Tv}(n) 
    # Fvj[i] = F[i,j] if i in find(F[j:n,j]), garbage otherwise
    
    nextp = Fp+1
    
    for j = 1:n
        for i in @view Fi[Fp[j]:Fp[j+1]-1]
            Fvj[i] = zero(Tv)
        end
        for p in Ap[j]:Ap[j+1]-1
            Fvj[Ai[p]] = Av[p]
        end
        @inbounds for k = @view Fk[Fq[j]:Fq[j+1]-1]
            Fvkj = Fv[Fp[k]]*Fv[nextp[k]]  # Factor out for performance
            for p = nextp[k]:Fp[k+1]-1
                # We compute a few dropped fill-ins here. It turns out computing 
                # and discarding is faster than introducing a branch. 
                Fvj[Fi[p]] -= Fv[p]*Fvkj
            end
            nextp[k] += 1
        end
        Fv[Fp[j]] = Fvj[j]
        for p = Fp[j]+1:Fp[j+1]-1
            Fv[p] = Fvj[Fi[p]]/Fvj[j]
        end
    end
    return Fv
end

end # module
