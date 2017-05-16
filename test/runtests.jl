using IncompleteSelectedInversion
using Base.Test

T = IncompleteSelectedInversion

@testset "iterate_jkp" begin
    srand(42)
    for i = 1:100
        n = rand(1:100)
        fill = rand(1:20)
        A = I + sprand(n,n,min(1.,fill/n))
        Ap,Ai = A.colptr,A.rowval

        jvals = Int[]
        iter = T.iterate_jkp(Ap,Ai)
        for (j,kpvals) in iter
            push!(jvals,j)
            kvals = Int[]
            for (k,pvals) in kpvals
                push!(kvals,k)
                ivals = Int[]
                for p in pvals
                    push!(ivals,Ai[p])
                end
                @test ivals == j-1+find(A[j:end,k])
            end
            @test sort(kvals) == find(A[j,1:j-1])
        end
        @test sort(jvals) == collect(1:n)
    end
end

@testset "symbolic" begin
    srand(42)
    for i = 1:100
        n = rand(1:100)
        fill = rand(1:20)
        A = I + sprand(n,n,min(1.,0.5*fill/n)); A += A'
        Ap,Ai = A.colptr,A.rowval

        Fp,Fi = T.symbolic(Ap,Ai,n)
        F = SparseMatrixCSC(n,n,Fp,Fi,ones(Bool,length(Fi)))
        F̂ = tril(lufact(full(A),Val{false}).factors .!= 0)
        @test (F == F̂) == true
    end
end

