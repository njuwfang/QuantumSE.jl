using QuantumSE
using Z3

ctx = Z3.Context()

@qprog id_circuit (n) begin
    for j in 1:n
        H(j)
    end
    for j in 1:n
        X(j)
    end
    for j in 1:n
        X(j)
    end
    for j in 1:n
        H(j)
    end
end

function prove_id(n)
    stabilizer = zeros(Bool, n, 2*n)
	phases = Vector{Z3.ExprAllocated}(undef, n)

    for j in 1:n
        stabilizer[j,j+n] = true
        phases[j] = _bv_const(ctx, "z_$(j)")
    end

    ρ01 = from_stabilizer(n, stabilizer, phases, ctx)
    ρ1 = copy(ρ01)

    stabilizer = zeros(Bool, n, 2*n)
	phases = Vector{Z3.ExprAllocated}(undef, n)

    for j in 1:n
        stabilizer[j,j] = true
        phases[j] = _bv_const(ctx, "x_$(j)")
    end

    ρ02 = from_stabilizer(n, stabilizer, phases, ctx)
    ρ2 = copy(ρ02)

    σ = CState([(:n, n),
            (:id_circuit, id_circuit),
            (:ctx, ctx),
        ])

    cfg1 = SymConfig(id_circuit(n), σ, ρ1)
    cfg2 = SymConfig(id_circuit(n), σ, ρ2)

    cfgs1 = QuantSymEx(cfg1)
    cfgs2 = QuantSymEx(cfg2)

    res = true
    for cfg in cfgs1
        if !check_state_equivalence(
            cfg.ρ, ρ01, (bool_val(ctx,true), cfg.ϕ[1], cfg.ϕ[2]);
            use_z3=true
           )
            res = false
            break
        end
    end

    if res
        for cfg in cfgs2
            if !check_state_equivalence(
                cfg.ρ, ρ02, (bool_val(ctx,true), cfg.ϕ[1], cfg.ϕ[2]);
                use_z3=true
               )
                res = false
                break
            end
        end
    end

    println(res)
end

prove_id(1)

for j in 1:20
    start = time()
    prove_id(j)
    println("Runtime for $(j): ", time() - start)
end
