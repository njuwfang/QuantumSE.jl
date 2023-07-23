using QuantumSE
using Z3
using SparseArrays
using LinearAlgebra: nullspace

ctx = Z3.Context()

struct PrimeG <: AbstractGroup
    p::Integer
    v::Integer

    PrimeG(p, v) = new(p, mod(v, p))
end

Base.:*(x::PrimeG, y::PrimeG) = PrimeG(x.p, x.v + y.v)
Base.inv(x::PrimeG) = PrimeG(x.p, x.p - x.v)

Hamming743 = GF2.([1 1 0 1 1 0 0;1 0 1 1 0 1 0;0 1 1 1 0 0 1])
Hamming733 = GF2.([1 0 0 0 1 1 0;0 1 0 0 1 0 1;0 0 1 0 0 1 1;0 0 0 1 1 1 1])

function TannerCode(G::Vector{<:AbstractGroup}, A::Vector{<:AbstractGroup}, B::Vector{<:AbstractGroup}, HA, HAt, HB, HBt)
    ng = length(G)
    na = length(A)
    nb = length(B)

    #  Cayley complex 
    F = [((g,0,0), (a*g,0,1), (a*g*b,1,1), (g*b,1,0)) for g in G for a in A for b in B]

    nf = length(F)

    # Tanner code
    F_idxs = Dict([(F[j], j) for j in 1:nf])

    IaX = Vector{Int}(undef, na*nb*ng*2)
    IbX = Vector{Int}(undef, na*nb*ng*2)
    IgX = Vector{Int}(undef, na*nb*ng*2)
    IfX = Vector{Int}(undef, na*nb*ng*2)
    IaZ = Vector{Int}(undef, na*nb*ng*2)
    IbZ = Vector{Int}(undef, na*nb*ng*2)
    IgZ = Vector{Int}(undef, na*nb*ng*2)
    IfZ = Vector{Int}(undef, na*nb*ng*2)
    @inbounds @simd for jg in 1:ng
        g = G[jg]
        for jb in 1:nb
            b = B[jb]
            for ja in 1:na
                a = A[ja]
                j = (jg-1)*na*nb*2+(jb-1)*na*2+(ja-1)*2+1
                IaX[j], IaX[j+1], IaZ[j], IaZ[j+1] = ja, ja, ja, ja
                IbX[j], IbX[j+1], IbZ[j], IbZ[j+1] = jb, jb, jb, jb
                IgX[j], IgX[j+1], IgZ[j], IgZ[j+1] = jg, jg+ng, jg, jg+ng
                IfX[j], IfX[j+1], IfZ[j], IfZ[j+1] = F_idxs[((g,0,0), (a*g,0,1), (a*g*b,1,1), (g*b,1,0))], F_idxs[((inv(a)*g*inv(b),0,0), (g*inv(b),0,1), (g,1,1), (inv(a)*g,1,0))], F_idxs[((inv(a)*g,0,0), (g,0,1), (g*b,1,1), (inv(a)*g*b,1,0))], F_idxs[((g*inv(b),0,0), (a*g*inv(b),0,1), (a*g,1,1), (g,1,0))]
            end
        end
    end

    IX = @. (IgX-1)*nf + IfX
    JX = @. (IaX-1)*nb + IbX
    VX = ones(GF2, na*nb*ng*2)

    IZ = @. (IgZ-1)*nf + IfZ
    JZ = @. (IaZ-1)*nb + IbZ
    VZ = ones(GF2, na*nb*ng*2)

    HXt = sparse(IX, JX, VX, nf*2*ng, na*nb)
    HZt = sparse(IZ, JZ, VZ, nf*2*ng, na*nb)

    HXt = reshape(HXt * kron(sparse(transpose(HA)), sparse(transpose(HB))), nf, :)
    HZt = reshape(HZt * kron(sparse(transpose(HAt)), sparse(transpose(HBt))), nf, :)

    HXt, HZt
end

function css_check(d, s, s_type, nq, adj)

    ## pre-condition
    ϕ₁ = bool_val(ctx, true)

    ## post-condition
    ϕ₂ = bool_val(ctx, true)
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:nq]
    for j in 1:length(s)
        ϕ₂ = ϕ₂ & (s[j] ⊻ reduce(⊻,  r[adj(j)]) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(nq)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(nq)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

@qprog tanner_x_m (idx) begin
    b = _xadj(idx)

    nb = length(b)
    for j in 2:nb
        CNOT(b[1], b[j])
    end
    H(b[1])
    res = M(b[1])
    H(b[1])
    for j in nb:-1:2
        CNOT(b[1], b[j])
    end

    res
end

@qprog tanner_z_m (idx) begin
    b = _zadj(idx)

    nb = length(b)
    for j in 2:nb
        CNOT(b[j], b[1])
    end
    res = M(b[1])
    for j in nb:-1:2
        CNOT(b[j], b[1])
    end

    res
end

@qprog tanner_decoder (nx, nz, nq, d) begin
    s_x = [tanner_x_m(j) for j in 1:nx]
    s_z = [tanner_z_m(j) for j in 1:nz]

    r_x = css_check(d, s_x, "X", nq, _xadj)
    r_z = css_check(d, s_z, "Z", nq, _zadj)

    for j in 1:nq
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    e = reduce(&, r_z[1:((d-1)÷2)])

    sX(1, e)
end

function check_tanner_decoder(m,k)
    #=
    p = 13
    q = 5

    G = pgl2(q)
    S = jacobi4squares(p, q)
    A = S
    B = S

    HA = [Hamming743 Hamming743]
    HAt = [Hamming733 Hamming733]
    HB = [Hamming733 Hamming733]
    HBt = [Hamming743 Hamming743]
    =#
    @info "Initailization Stage"
    @time begin

        #G = [PrimeG(7^m*2^k, j-1) for j in 1:7^m*2^k]
        #A = [PrimeG(7^m*2^k, (j-1)*7^(m-1)*2^k) for j in 1:7]
        #B = [PrimeG(7^m*2^k, (j-1)*7^(m-1)*2^k) for j in 1:7]
        aa = k
        G = [PrimeG(7^m*aa, j-1) for j in 1:7^m*aa]
        A = [PrimeG(7^m*aa, (j-1)*7^(m-1)*aa) for j in 1:7]
        B = [PrimeG(7^m*aa, (j-1)*7^(m-1)*aa) for j in 1:7]

        #HA = Hamming313
        #HAt = Hamming322
        #HB = Hamming322
        #HBt = Hamming313
        
        HA = Hamming743
        HAt = Hamming733
        HB = Hamming733
        HBt = Hamming743

        HXt, HZt = TannerCode(G, A, B, HA, HAt, HB, HBt)

        n, nx = size(HXt)
        nz = size(HZt,2)

        @show nx, nz

        X_idxs = [findall(!iszero, HXt[:,j]) for j in 1:nx]
        Z_idxs = [findall(!iszero, HZt[:,j]) for j in 1:nz]

        _xadj(i) = X_idxs[i]
        _zadj(i) = Z_idxs[i]

        ρ01, ρ02, dx, dz = from_css_code(Matrix{GF2}(transpose(HXt)), Matrix{GF2}(transpose(HZt)), ctx)
        d = 6 #min(dx, dz)

        ρ1 = copy(ρ01)
        ρ2 = copy(ρ02)

        σ = CState([
            (:d, d),
            (:tanner_decoder, tanner_decoder),
            (:tanner_x_m, tanner_x_m),
            (:tanner_z_m, tanner_z_m),
            (:_xadj, _xadj),
            (:_zadj, _zadj),
            (:ctx, ctx),
            (:css_check, css_check)
        ])

        num_x_errors = (d-1)÷2
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, n) == bv_val(ctx, num_x_errors, _len2(n)+1)
        
        #num_z_errors = (d-1)÷2
        #z_errors = inject_errors(ρ1, "Z")
        #ϕ_z1 = _sum(ctx, z_errors, n) == bv_val(ctx, num_z_errors, _len2(n)+1)
        
        x_errors = inject_errors(ρ2, "X")
        ϕ_x2 = _sum(ctx, x_errors, n) == bv_val(ctx, num_x_errors, _len2(n)+1)

        #num_z_errors = (d-1)÷2
        #z_errors = inject_errors(ρ2, "Z")
        #ϕ_z2 = _sum(ctx, z_errors, n) == bv_val(ctx, num_z_errors, _len2(n)+1)

        cfg1 = SymConfig(tanner_decoder(nx,nz,n,d), σ, ρ1)
        cfg2 = SymConfig(tanner_decoder(nx,nz,n,d), σ, ρ2)
    end

    @info "Symbolic Execution Stage"
    @time begin
        cfgs1 = QuantSymEx(cfg1)
        cfgs2 = QuantSymEx(cfg2)
    end

    @info "SMT Solver Stage"
    @time begin
        res = true
        for cfg in cfgs1
            if !check_state_equivalence(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg.ϕ[1], cfg.ϕ[2]),
                `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
                res = false
                break
            end
        end

        if res
            for cfg in cfgs2
                if !check_state_equivalence(
                    cfg.ρ, ρ02, (ϕ_x2 #=& ϕ_z2=#, cfg.ϕ[1], cfg.ϕ[2]),
                    `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
                    res = false
                    break
                end
            end
        end
    end

    res
end
