using QSE
using Z3

ctx = Context()

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

    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    adj = s_type == "X" ? _xadj : _zadj
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:2*d*d]
    for j in 1:d*d
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[adj(d, j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

function toric_x_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in _xadj(d, idx)
		s[j] = true
	end

	s
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

function toric_z_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in (_zadj(d, idx) .+ 2*d*d)
		s[j] = true
	end

	s
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

function toric_lx1(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[i*d] = true
	end

	s
end

function toric_lx2(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[d*d+i] = true
	end

	s
end

function toric_lz1(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[3*d*d+i*d] = true
	end

	s
end

function toric_lz2(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[2*d*d+i] = true
	end

	s
end

@qprog toric_decoder (d) begin

    s_x = [toric_x_m(d, j) for j in 1:d*d]
    s_z = [toric_z_m(d, j) for j in 1:d*d]

    s_z[d*d*2÷3] = _bv_val(ctx, 0) # a strange bug
    
    r_x = mwpm(d, s_x, "X")
    r_z = mwpm(d, s_z, "Z")
    
    for j in 1:2*d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

end

function check_toric_decoder(d::Integer)

	num_qubits = d*d*2

	stabilizer = Matrix{Bool}(undef, num_qubits, 2*num_qubits)
	phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
	lx = _bv_const(ctx, "lx")
	lz = _bv_const(ctx, "lz")

	@simd for i in 1:d*d-1
		stabilizer[i,:] = toric_x_s(d, i)
		stabilizer[i+d*d,:] = toric_z_s(d, i)
		phases[i] = _bv_val(ctx, 0)
		phases[i+d*d] = _bv_val(ctx, 0)
	end

	stabilizer[d*d,:] = toric_lx1(d)
	phases[d*d] = lx
	stabilizer[2*d*d,:] = toric_lz1(d)
	phases[2*d*d] = lz

    σ = CState([(:d, d),
        (:toric_code, toric_decoder),
        (:toric_z_m, toric_z_m),
        (:toric_x_m, toric_x_m),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:mwpm, mwpm)
    ])

    ρ₀ = from_stabilizer(num_qubits, stabilizer, phases, ctx)
    ρ = QState(ρ₀)

    ϕ_x = inject_errors(ρ, (d-1)÷2, "X")
	ϕ_z = inject_errors(ρ, (d-1)÷2, "Z")

    cfg0 = SymConfig(toric_decoder(d), σ, ρ)

    for cfg in QuantSymEx(cfg0)
        check_state_equivalence(
            cfg.ρ, ρ₀,
            (ϕ_x & ϕ_z, cfg.ϕ[1], cfg.ϕ[2]),
            `yices-smt2`
        )
    end

end