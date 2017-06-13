using IncompleteSelectedInversion
using Base.Test

T = IncompleteSelectedInversion

@testset "errors" begin
    Ap = [1,2]
    Ai = [2]
    @test_throws Exception T.checkmat(Ap,Ai)

    Ap = [1,3]
    Ai = [1]
    @test_throws Exception T.checkmat(Ap,Ai)

    Ap = [1,3]
    Ai = [1,2]
    Av = [1]
    @test_throws Exception T.checkmat(Ap,Ai,Av)
end

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

@testset "numeric" begin
    srand(42)
    for i = 1:100
        n = rand(1:100)
        fill = rand(1:20)
        A = I + sprand(n,n,min(1.,0.5*fill/n)); A += A'
        Ap,Ai,Av = A.colptr,A.rowval,A.nzval

        Fp,Fi = T.symbolic(Ap,Ai,n)
        Fv = T.numeric(Ap,Ai,Av,Fp,Fi)
        F = SparseMatrixCSC(n,n,Fp,Fi,Fv)
        F̂ = tril(lufact(full(A),Val{false}).factors)
        @test (F ≈ F̂) == true
    end
end

@testset "selinv_jki" begin
    srand(42)
    for i = 1:100
        n = rand(1:100)
        fill = rand(1:20)
        A = 3I + sprand(n,n,min(1.,0.5*fill/n)); A += A'
        #   ^ need to make sure matrix is sufficiently well conditioned
        Ap,Ai,Av = A.colptr,A.rowval,A.nzval

        Fp,Fi = T.symbolic(Ap,Ai,n)
        Fv = T.numeric(Ap,Ai,Av,Fp,Fi)
        Bv = T.selinv_jki(Fp,Fi,Fv)
        B = SparseMatrixCSC(n,n,Fp,Fi,Bv)
        B̂ = inv(full(A))
        @test all(Bi == 0 || Bi ≈ B̂i for (Bi,B̂i) in zip(B,B̂))
    end
end

@testset "selinv_kij" begin
    srand(42)
    for i = 1:100
        n = rand(1:100)
        fill = rand(1:20)
        A = 3I + sprand(n,n,min(1.,0.5*fill/n)); A += A'
        #   ^ need to make sure matrix is sufficiently well conditioned
        Ap,Ai,Av = A.colptr,A.rowval,A.nzval

        Fp,Fi = T.symbolic(Ap,Ai,n)
        Fv = T.numeric(Ap,Ai,Av,Fp,Fi)
        Fq,Fj,Fw = T.permute4selinv(Fp,Fi,Fv)
        Bw = T.selinv_kij(Fq,Fj,Fw)
        B = T.packsparse(Fq,Fj,Bw)
        B̂ = inv(full(A))[end:-1:1,end:-1:1]
        @test all(Bi == 0 || Bi ≈ B̂i for (Bi,B̂i) in zip(B,B̂))
    end
end

