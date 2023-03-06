using QSE
using Z3
using Test

@testset "Toric_Decoder" begin

    function _xadj(d, idx)
        ai = (idx-1)÷d
        bi = (idx-1)%d
    
        [d*d+ai*d+bi, ai*d+(bi+d-1)%d, ai*d+bi, d*d+((ai+d-1)%d)*d+bi] .+ 1
    end
    
    function _zadj(d, idx)
        ai = (idx-1)÷d
        bi = (idx-1)%d
    
        [((ai+1)%d)*d+bi, d*d+ai*d+bi, d*d+ai*d+(bi+1)%d, ai*d+bi] .+ 1
    end

    function mwpm(d::Integer, s, s_type)

        ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

        ϕ₂ = bool_val(ctx, true)
        adj = s_type == "X" ? _xadj : _zadj
        r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:2*d*d]
        for j in 1:d*d
            ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[adj(d, j)...]])) == _bv_val(ctx, 0))
        end

        ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

        (r, ϕ₁ & ϕ₂ & ϕ₃)
    end
    
    @qprog toric_x_m (d, idx) begin
        b = _xadj(d, idx)
    
        CNOT(b[1], b[2])
        CNOT(b[3], b[4])
        CNOT(b[3], b[1])
        H(b[3])
        res = M(b[3])
        H(b[3])
        CNOT(b[3], b[1])
        CNOT(b[3], b[4])
        CNOT(b[1], b[2])
    
        res
    end
    
    @qprog toric_z_m (d, idx) begin
        b = _zadj(d, idx)
    
        CNOT(b[1], b[2])
        CNOT(b[3], b[4])
        CNOT(b[2], b[4])
        res = M(b[4])
        CNOT(b[2], b[4])
        CNOT(b[3], b[4])
        CNOT(b[1], b[2])
    
        res
    end
    
    @qprog toric_decoder (d) begin
        s_x = [toric_x_m(d, j) for j in 1:d*d]
        s_z = [toric_z_m(d, j) for j in 1:d*d]
    
        r_x = mwpm(d, s_x, "X")
        r_z = mwpm(d, s_z, "Z")
    
        for j in 1:2*d*d
            sZ(j, r_x[j])
            sX(j, r_z[j])
        end
    end

    @qprog test () begin
        x = [1, 2, 3]
        x[2] = 100
    end

    ctx = Z3.Context()

    d = 5
    σ = CState([(:d, d),
        (:toric_code, toric_decoder),
        (:toric_z_m, toric_z_m),
        (:toric_x_m, toric_x_m),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:mwpm, mwpm)
    ])

    ρ = QState(2*d*d, ctx)

    cfg = SymConfig(test(), σ, ρ)

    println(QuantSymEx(cfg)[1].σ)


end