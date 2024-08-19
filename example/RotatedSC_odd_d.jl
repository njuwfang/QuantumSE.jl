using QuantumSE
using Z3

N4 = 4
N2 = 2

AD4 = 4
AD2_H = 2
AD2_V = 2

ctx = Context()

function _xadj(d, idx, nbr)
    r = (idx-1)÷d
    c = (idx-1)%d

    if nbr == 4
        return [r*d+c, r*d+c+1, (r+1)*d+c, (r+1)*d+(c+1)] .+ 1
    elseif nbr == 2
        return [r*d + c, r*d + c+1 ] .+ 1
    else
        return [] 
    end

end

function _zadj(d, idx, nbr)
    r = (idx-1)÷d
    c = (idx-1)%d

    if nbr == 4
        return [r*d+c, r*d+c+1, (r+1)*d+c, (r+1)*d+(c+1)] .+ 1
    elseif nbr == 2
        return [r*d + c, (r+1)*d + c] .+ 1
    else
        return [] 
    end

end


function mwpm(d::Integer, s, s_type, x_4q, x_2q, z_4q, z_2q)


    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)


    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:d*d]

    if s_type == "X"
        for j in eachindex(x_4q)
            ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_xadj(d, x_4q[j], 4)...]])) == _bv_val(ctx, 0))
        end
        for j in eachindex(x_2q)
            ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_xadj(d, x_2q[j], 2)...]])) == _bv_val(ctx, 0))
        end
    else
        for j in eachindex(z_4q)
            ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_xadj(d, z_4q[j], 4)...]])) == _bv_val(ctx, 0))
        end
        for j in eachindex(z_2q)
            ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_xadj(d, z_2q[j], 2)...]])) == _bv_val(ctx, 0))
        end
    end


    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

@qprog surface_code_z_m (d, idx, nbr) begin
    b = _zadj(d, idx, nbr)

    if length(b) == 4
        CNOT(b[1], b[2])
        CNOT(b[3], b[4])
        CNOT(b[2], b[4])
        res = M(b[4])
        CNOT(b[2], b[4])
        CNOT(b[3], b[4])
        CNOT(b[1], b[2])
    else
        if length(b) == 2
            CNOT(b[1], b[2])
            res = M(b[1])
            CNOT(b[1], b[2])
        else
            res = -100
        end
    end

    res
end

@qprog surface_code_x_m (d, idx, nbr) begin
    b = _xadj(d, idx, nbr)

    if length(b) == 4
        CNOT(b[1], b[2])
        CNOT(b[3], b[4])
        CNOT(b[3], b[1])
        H(b[3])
        res = M(b[3])
        H(b[3])
        CNOT(b[3], b[1])
        CNOT(b[3], b[4])
        CNOT(b[1], b[2])

    else
        if length(b) == 2
            CNOT(b[1], b[2])
            H(b[1])
            res = M(b[1])
            H(b[1])
            CNOT(b[1], b[2])
        else
            res = -100
        end
    end

    res
end


@qprog surface_code_decoder (d, x_4q, x_2q, z_4q, z_2q) begin
    # print("Decoder start")

    s_z4 = [surface_code_z_m(d, idx, 4) for idx in z_4q]
    s_x4 = [surface_code_x_m(d, idx, 4) for idx in x_4q]

    s_z2 = [surface_code_z_m(d, idx, 2) for idx in z_2q]
    s_x2 = [surface_code_x_m(d, idx, 2) for idx in x_2q]

    r_x = mwpm(d, vcat(s_x2, s_x4), "X", x_4q, x_2q, z_4q, z_2q)
    r_z = mwpm(d, vcat(s_z2, s_z4), "Z", x_4q, x_2q, z_4q, z_2q)

    # print("rx=$(r_x)")
    # print("rz=$(r_z)")


    for j in 1:d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    # e = reduce(&, r_z[1:(d-1)÷2])

    # sX(1, e)

    # print("Decoder end")
end



function check_surface_code_decoder(d::Integer)

    d = 3

    @info "Initialization Stage"
    t0 = time()
    begin
        num_qubits = d*d

	    stabilizer = fill(false, num_qubits, 2*num_qubits)

        z_4q = [ 1 + r*d+c+r%2 for r in 0:(d-2) for c in 0:2:(d-2) ]
        # println("Z 4q: $(z_4q)")

        x_4q = [ 1 + r*d+c+(r+1)%2 for r in 0:(d-2) for c in 0:2:(d-2) ]
        # println("X 4q: $(x_4q)")


        x_2q_up = [ 1+c for c in 0:2:(d-2) ]
        x_2q_down = [ 1+ (d-1)*d + c for c in 1:2:(d-2) ]
        x_2q = vcat(x_2q_up, x_2q_down)

    #    println("X 2q: $(x_2q)")


        z_2q_right = [ 1 + r*d+ (d-1) for r in 0:2:(d-2) ]
        z_2q_left = [ 1 + r*d for r in 1:2:(d-2) ]
        z_2q = vcat(z_2q_left, z_2q_right)

    #    println("Z 2q: $(z_2q)")



        r = 1

        for j in eachindex(x_4q)
            for x_adj in  _xadj( d, x_4q[j], 4)
                stabilizer[r, x_adj] = true
            end
            r += 1
        end



        for j in eachindex(x_2q)
            for x_adj in _xadj( d, x_2q[j], 2)
                stabilizer[r, x_adj] = true
            end
            
            r += 1
        end


        for j in eachindex(z_4q)
            for z_adj in _zadj( d, z_4q[j], 4)
                stabilizer[r, num_qubits+z_adj] = true
            end
            r += 1
        end



        for j in eachindex(z_2q)
            for z_adj in _zadj( d, z_2q[j], 2)
                stabilizer[r, num_qubits+z_adj] = true
            end
            r += 1
        end

        # r: 1, 2, ..., d ;; c = 1
        sc_lx = [ 1+ r*d + 1 for r in 0:(d-1) ]

        for idx in sc_lx
            stabilizer[num_qubits, idx] = true
        end
        
        # println("Encoded stabilizer : $(stabilizer)")



	    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
	    lx = _bv_const(ctx, "lx")
	    lz = _bv_const(ctx, "lz")

	    for i in 1:num_qubits-1
	    	phases[i] = _bv_val(ctx, 0)
	    end

        phases[num_qubits] = lx


        ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ1 = copy(ρ01)


        σ = CState([(:d, d),
            (:surface_code_decoder, surface_code_decoder),
            (:surface_code_z_m, surface_code_z_m),
            (:surface_code_x_m, surface_code_x_m),
            (:_xadj, _xadj),
            (:_zadj, _zadj),
            (:ctx, ctx),
            (:mwpm, mwpm)
        ])

        num_x_errors = (d-1)÷2
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

	  
        decoder = surface_code_decoder(d, x_4q, x_2q, z_4q, z_2q)
        cfg1 = SymConfig(decoder, σ, ρ1)
        # cfg2 = SymConfig(surface_code_decoder(d), σ, ρ2)
    end

    @info "Symbolic Execution Stage"
    t1 = time()
    begin
        cfgs1 = QuantSymEx(cfg1)
        # cfgs2 = QuantSymEx(cfg2)
    end

    @info "SMT Solver Stage"
    t2 = time()

    begin
        res = true
        # for cfg in cfgs1
        if !check_state_equivalence(
            ρ01, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg1.ϕ[1], cfg1.ϕ[2]),
            `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`
            #`bitwuzla --smt-comp-mode true -S kissat`
            #`bitwuzla --smt-comp-mode true -rwl 0`
            )
            res = false
            # break
        end
        # end

    end

    t3 = time()

    res, t2-t0, t1-t0, t2-t1, t3-t2
end

# check_surface_code_decoder(3) # precompile time
open("surface_code.dat", "w") do io
    println(io, "d: res nq all init qse smt")
    for d in 3:2:15
        res_d, all, init, qse, smt = check_surface_code_decoder(d)
        println("d: res nq all init qse smt")
        println("$(d): $(res_d) $(d*d) $(all) $(init) $(qse) $(smt)")
        println(io, "$(d): $(res_d) $(d*d) $(all) $(init) $(qse) $(smt)")
    end
end
