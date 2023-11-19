using QuantumSE
using Z3

ctx = Context()

_adj(idx) = [idx, idx+1]

function mwpm(n::Integer, s)

    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    r = [_bv_const(ctx, "r_$(j)") for j in 1:n]
    for j in 1:n-1
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_adj(j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(n)), x)).(r) ) <= bv_val(ctx, (n-1)÷2, _len2(n)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

@qprog repetition_m_zz (n, idx) begin
    CNOT(idx, idx%n+1)
    res = M(idx%n+1)
    CNOT(idx, idx%n+1)
    
    res
end

function repetition_s(n, idx)
    s = zeros(Bool, 2*n)

    s[n+idx] = true
    s[n+idx+1] = true

    s
end

function repetition_lx(n)
    s = zeros(Bool, 2*n)

    s[n+1] = true

    s
end

function repetition_lz(n)
    s = zeros(Bool, 2*n)

    for j in 1:n
        s[j] = true
    end

    s
end

@qprog repetition_decoder (n) begin
    s = [repetition_m_zz(n, j) for j in 1:n]

    r = mwpm(n, s)

    for j in 1:n
        sX(j, r[j])
    end

    e = reduce(&, r[1:((n-1)÷2)])

    sX(1, e)
end

function check_repetition_decoder(n)
    @info "Initailization Stage"
    t0 = time()
    @time begin
        num_qubits = n

	    stabilizer = Matrix{Bool}(undef, num_qubits, 2*num_qubits)
	    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
	    lx = _bv_const(ctx, "lx")
	    lz = _bv_const(ctx, "lz")

	    @simd for i in 1:n-1
	    	stabilizer[i+1,:] = repetition_s(n, i)
	    	phases[i+1] = _bv_val(ctx, 0)
	    end

	    stabilizer[1,:] = repetition_lx(n)
	    phases[1] = lx

        σ = CState([(:n, n),
            (:repetition_decoder, repetition_decoder),
            (:repetition_m_zz, repetition_m_zz),
            (:_adj, _adj),
            (:ctx, ctx),
            (:mwpm, mwpm)
        ])

        ρ₀ = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ = copy(ρ₀)

        num_x_errors = (n-1)÷2
        x_errors = inject_errors(ρ, "X")
        ϕ_x = _sum(ctx, x_errors, num_qubits) <= bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

	    #num_z_errors = 0
        #z_errors = inject_errors(ρ, "Z")
        #ϕ_z = _sum(ctx, z_errors, num_qubits) <= bv_val(ctx, num_z_errors, _len2(num_qubits)+1)

        cfg0 = SymConfig(repetition_decoder(n), σ, ρ)
    end

    @info "Symbolic Execution Stage"
    t1 = time()
    @time cfgs = QuantSymEx(cfg0)

    @info "SMT Solver Stage"
    t2 = time()
    @time begin
        res = true
        for cfg in cfgs
            if !check_state_equivalence(
                cfg.ρ, ρ₀, (ϕ_x #=& ϕ_z=#, cfg.ϕ[1], cfg.ϕ[2]),
                `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
                res = false
                break
            end
        end
    end

    t3 = time()

    res, t3-t0, t1-t0, t2-t1, t3-t2
end

check_repetition_decoder(20) # precompile time

open("repetition_code.dat", "w") do io
  println(io, "nq all init qse smt")
  println("nq all init qse smt")
  for j in 1:28
    res, all, init, qse, smt = check_repetition_decoder(50*j)
    println(io, "$(50*j) $(all) $(init) $(qse) $(smt)")
    println("$(j)/28: $(50*j) $(all) $(init) $(qse) $(smt)")
  end
end
