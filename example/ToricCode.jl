using QuantumSE
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

    println("precondition 1: $(ϕ₁)")

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    adj = s_type == "X" ? _xadj : _zadj
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:2*d*d]
    for j in 1:d*d
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[adj(d, j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

    println("post-condition 2: $(ϕ₂)")
    println("post-condition 3: $(ϕ₃)")

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

    # d, 2d, 3d, ... , d*d
	@inbounds @simd for i in 1:d
		s[i*d] = true
	end

	s
end

function toric_lx2(d::Integer)
	s = zeros(Bool, 4*d*d)

    # d*d + [1, 2, ..., d]
	@inbounds @simd for i in 1:d
		s[d*d+i] = true
	end

	s
end

function toric_lz1(d::Integer)
	s = zeros(Bool, 4*d*d)

    # Z-> d*d + [d, 2d, 3d, ... , d*d]
	@inbounds @simd for i in 1:d
		s[2*d*d + d*d+i*d] = true
	end

	s
end

function toric_lz2(d::Integer)
	s = zeros(Bool, 4*d*d)

    # Z-> [1, 2, ..., d]
	@inbounds @simd for i in 1:d
		s[2*d*d+i] = true
	end

	s
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

    # a strange bug
    e = reduce(&, r_z[1:(d-1)÷2])

    sX(1, e)

end

function check_toric_decoder(d::Integer)

    @info "Initailization Stage"
    t0 = time()
    begin
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

        ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ1 = copy(ρ01)

        stabilizer[d*d,:] = toric_lx2(d)
	    phases[d*d] = lx
	    stabilizer[2*d*d,:] = toric_lz2(d)
	    phases[2*d*d] = lz

        ρ02 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ2 = copy(ρ02)

        σ = CState([(:d, d),
            (:toric_decoder, toric_decoder),
            (:toric_z_m, toric_z_m),
            (:toric_x_m, toric_x_m),
            (:_xadj, _xadj),
            (:_zadj, _zadj),
            (:ctx, ctx),
            (:mwpm, mwpm)
        ])

        num_x_errors = (d-1)÷2
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

	    #num_z_errors = (d-1)÷2
        #z_errors = inject_errors(ρ1, "Z")
        #ϕ_z1 = _sum(ctx, z_errors, num_qubits) == bv_val(ctx, num_z_errors, _len2(num_qubits)+1)

        x_errors = inject_errors(ρ2, "X")
        ϕ_x2 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

	    #num_z_errors = (d-1)÷2
        #z_errors = inject_errors(ρ2, "Z")
        #ϕ_z2 = _sum(ctx, z_errors, num_qubits) == bv_val(ctx, num_z_errors, _len2(num_qubits)+1)

        cfg1 = SymConfig(toric_decoder(d), σ, ρ1)
        cfg2 = SymConfig(toric_decoder(d), σ, ρ2)
    end

    @info "Symbolic Execution Stage"
    t1 = time()
    begin
        cfgs1 = QuantSymEx(cfg1)
        cfgs2 = QuantSymEx(cfg2)
    end

    @info "SMT Solver Stage"
    t2 = time()
    begin
        res = true
        for cfg in cfgs1
            if !check_state_equivalence(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg.ϕ[1], cfg.ϕ[2]),
                `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`
                #`bitwuzla --smt-comp-mode true -S kissat`
                #`bitwuzla --smt-comp-mode true -rwl 0`
               )
                res = false
                break
            end
        end

        if res
            for cfg in cfgs2
                if !check_state_equivalence(
                    cfg.ρ, ρ02, (ϕ_x2 #=& ϕ_z2=#, cfg.ϕ[1], cfg.ϕ[2]),
                    `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`
                    #`bitwuzla --smt-comp-mode true -S kissat`
                    #`bitwuzla --smt-comp-mode true -rwl 0`
                   )
                    res = false
                    break
                end
            end
        end
    end

    t3 = time()

    res, t3-t0, t1-t0, t2-t1, t3-t2
end

# check_toric_decoder(3) # precompile time

open("toric_code.dat", "w") do io
  println(io, "nq res all init qse smt")
  println("nq res all init qse smt")
  for d in 3:3
    res, all, init, qse, smt = check_toric_decoder(d)
    println(io, "$(d) : $(res) $(d*d*2) $(all) $(init) $(qse) $(smt)")
    println("$(d): $(res) $(d*d*2) $(all) $(init) $(qse) $(smt)")
  end
end
