module QuantumSE

abstract type AbstractSymQuantumState end

include("QuantumProgram.jl")
export @qprog, QProg, QEmpty

include("GF2.jl")
export GF2

include("SymbolicStabilizer.jl")
export M, X, Y, Z, sX, sY, sZ, S, H, CNOT,
    SymStabilizerState, from_stabilizer, print_full_tableau, update!, inject_errors, from_css_code,
    _bv_val, _bv_const, _len2, check_state_equivalence, _sum

include("SymbolicExecution.jl")
export CState, SymConfig, QuantSymEx, CEval

include("LinearGroup.jl")
export AbstractGroup, AbstraLinearGroup, L2, PL2, pgl2, psl2, jacobi4squares

include("Sampler.jl")
export Sampler, parse_stim

using PrecompileTools
using Z3

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    _ctx = Z3.Context()
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
    n = 4
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        stabilizer = zeros(Bool, n, 2*n)
	    phases = Vector{Z3.ExprAllocated}(undef, n)

        for j in 1:n
            stabilizer[j,j+n] = true
            phases[j] = _bv_const(_ctx, "z_$(j)")
        end

        ρ0 = from_stabilizer(n, stabilizer, phases, _ctx)
        ρ = copy(ρ0)

        σ = CState([(:n, n),
                (:id_circuit, id_circuit),
                (:ctx, _ctx),
            ])

        cfg = SymConfig(id_circuit(n), σ, ρ)

        cfgs = QuantSymEx(cfg)

        res = true
        for cfg in cfgs
            if !check_state_equivalence(
                cfg.ρ, ρ0, (bool_val(_ctx,true), cfg.ϕ[1], cfg.ϕ[2]);
                use_z3=true
                )
                res = false
                break
            end
        end
    end
end

end # module QSE
