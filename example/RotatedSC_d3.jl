using QuantumSE
using Z3

ctx = Context()

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


@qprog surface_code_x_m (idx) begin
    b = _xadj(idx)

    println("SC X_m($(idx))...")

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

@qprog surface_code_z_m (idx) begin
    b = _zadj(idx)

    println("SC Z_m($(idx))...")


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

@qprog surface_code_decoder (d) begin

    nq = d*d 

    s_x = [surface_code_x_m(j) for j in 1:4]
    s_z = [surface_code_z_m(j) for j in 1:4]

    r_x = css_check(d, s_x, "X", nq, _xadj)
    r_z = css_check(d, s_z, "Z", nq, _zadj)

    # print("rx=$(r_x)")
    # print("rz=$(r_z)")


    for j in 1:d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    e = reduce(&, r_z[1:(d-1)÷2])

    sX(1, e)

end



function check_surface_code_decoder(d::Integer)

    d = 3

    @info "Initailization Stage"
    t0 = time()
    begin
        num_qubits = d*d

	    stabilizer = fill(false, num_qubits, 2*num_qubits)

        X_idxs = [[2, 3, 5, 6], [4, 5, 7, 8], [1, 2], [8, 9]]

        Z_idxs = [[1, 2, 4, 5], [5, 6, 8, 9], [4, 7], [3, 6]]


        #_xadj(i) = [ x for x in X_idxs if x[1] == i][1] # _xadj(2) = Vector(2, 3, 5, 6)
        _xadj(j) = X_idxs[j]
        #_zadj(i) = [ z for z in Z_idxs if z[1] == i][1]
        _zadj(j) = Z_idxs[j]

        #xq = [1, 2, 4, 8]

        for r in 1:4
            for x_adj in X_idxs[r]#_xadj(xq[r])
                # println("r=$(r), x_adj=$(x_adj)")
                stabilizer[r, x_adj] = true
            end
        end

        #zq = [1, 3, 4, 5]

        for r in 1:4
            for z_adj in Z_idxs[r]#_zadj(zq[r])
                stabilizer[4+r, num_qubits+z_adj] = true
            end
        end

        stabilizer[9,2] = true
        stabilizer[9,5] = true
        stabilizer[9,8] = true

        
        
        println("Encoded stabilizer : $(stabilizer)")



	    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
	    lx = _bv_const(ctx, "lx")
	    lz = _bv_const(ctx, "lz")

	    for i in 1:num_qubits-1
	    	phases[i] = _bv_val(ctx, 0)
	    end

        phases[9] = lx

        ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ1 = copy(ρ01)


        σ = CState([(:d, d),
            (:surface_code_decoder, surface_code_decoder),
            (:surface_code_z_m, surface_code_z_m),
            (:surface_code_x_m, surface_code_x_m),
            (:_xadj, _xadj),
            (:_zadj, _zadj),
            (:ctx, ctx),
            (:css_check, css_check)
        ])

        num_x_errors = (d-1)÷2
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

	  
        decoder = surface_code_decoder(d)
        cfg1 = SymConfig(decoder, σ, ρ1)
        # cfg2 = SymConfig(surface_code_decoder(d), σ, ρ2)
    end

    @info "Symbolic Execution Stage"
    t1 = time()
    begin
        cfgs1 = QuantSymEx(cfg1) # ERROR: LoadError: BoundsError: attempt to access 0-element Vector{Vector{Int64}} at index [1]
        # cfgs2 = QuantSymEx(cfg2)
    end

    @info "SMT Solver Stage"
    t2 = time()

    begin
        res = true
        for cfg in cfgs1
        if !check_state_equivalence(
            # ρ01, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg1.ϕ[1], cfg1.ϕ[2]),
            cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg.ϕ[1], cfg.ϕ[2]),
            `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`
            #`bitwuzla --smt-comp-mode true -S kissat`
            #`bitwuzla --smt-comp-mode true -rwl 0`
            )
            res = false
            break
        end
        end

    end

    t3 = time()

    res, t2-t0, t1-t0, t2-t1, t3-t2
end

# check_surface_code_decoder(3) # precompile time

d = 3
res, all, init, qse, smt = check_surface_code_decoder(d)
println("d: res nq all init qse smt")
println("$(d): $(res) $(d*d*2) $(all) $(init) $(qse) $(smt)")