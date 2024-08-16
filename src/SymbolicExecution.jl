using Z3
using MacroTools: postwalk, @capture

const CState = Dict{Symbol,Any}

const SymProbs = Dict{Z3.ExprAllocated,Z3.ExprAllocated}

mutable struct SymConfig
    S::QProg
    σ::CState
    const ρ::AbstractSymQuantumState
    P::SymProbs
    ϕ::Tuple{Z3.ExprAllocated, Z3.ExprAllocated}

    const ctx::Z3.ContextAllocated

    SymConfig(num_qubits::Integer) = begin
        ctx = Z3.Context()
        new(QEmpty, CState(), SymStabilizerState(num_qubits, ctx), SymProbs(), (bool_val(ctx, true), bool_val(ρ.ctx, true)), ctx, [])
    end

    SymConfig(S::QProg, σ::CState, ρ::AbstractSymQuantumState, P::SymProbs, ϕ::Tuple{Z3.ExprAllocated, Z3.ExprAllocated}, ctx::Z3.ContextAllocated) = new(copy(S), copy(σ), copy(ρ), copy(P), (ϕ[1], ϕ[2]), ctx)

    SymConfig(S::QProg, σ::CState, ρ::AbstractSymQuantumState) = SymConfig(S, σ, ρ, SymProbs(), (bool_val(ρ.ctx, true), bool_val(ρ.ctx, true)), ρ.ctx)

    SymConfig(cfg::SymConfig) = SymConfig(cfg.S, cfg.σ, cfg.ρ, cfg.P, (cfg.ϕ[1], cfg.ϕ[2]), cfg.ctx)
end

Base.copy(cfg::SymConfig) = SymConfig(cfg)

CEval(σ::CState, e) = eval(postwalk(x -> x isa Symbol ? x in keys(σ) ? σ[x] : x : x, e))

function CAssign(σ, lvalue, rvalue)
    if lvalue isa Expr
        if lvalue.head == :ref
            σ[lvalue.args[1]][CEval(σ, lvalue.args[2])] = CEval(σ, rvalue)
        end
    else
        σ[lvalue] = CEval(σ, rvalue)
    end
end

function QuantSymEx(cfg::SymConfig)

    length(cfg.S.args) != 0 || return [cfg]

    inst = cfg.S.args[1]
    if ~isa(inst, Expr)
        cfg.σ[:__res__] = CEval(cfg.σ, inst)
        cfg.S.args = cfg.S.args[2:end]
        return QuantSymEx(cfg)
    end

    if inst.head == :call
        
        if inst.args[1] in [:H, :S, :X, :Y, :Z, :CNOT] # Clifford gates
            if length(inst.args) == 2
                eval(inst.args[1])(cfg.ρ, CEval(cfg.σ, inst.args[2]))
            elseif length(inst.args) == 3
                eval(inst.args[1])(cfg.ρ, CEval(cfg.σ, inst.args[2]), CEval(cfg.σ, inst.args[3]))
            end
        elseif inst.args[1] in [:sX, :sY, :sZ] # Symbolic Pauli gates
            eval(inst.args[1])(cfg.ρ, CEval(cfg.σ, inst.args[2]), CEval(cfg.σ, inst.args[3]))
        elseif inst.args[1] == :M # Measurement without lvalue
            cfg.σ[:__res__] = eval(inst.args[1])(cfg.ρ, CEval(cfg.σ, inst.args[2]), eval(CEval(cfg.σ, inst.args[3])))
        else
            S = copy(cfg.S)
            res = CEval(cfg.σ, inst)
            if res isa Expr
                σ = cfg.σ
                cfg.S.args = CEval(cfg.σ, inst).args
                cfg = QuantSymEx(cfg)[1]
                σ[:__res__] = cfg.σ[:__res__]
                cfg.σ = σ
            elseif res isa Tuple
                cfg.σ[:__res__] = res[1]
                cfg.ϕ = (cfg.ϕ[1] & res[2], cfg.ϕ[2] & res[3])
            else
                cfg.σ[:__res__] = res
            end
            cfg.S.args = S.args[2:end]
            return QuantSymEx(cfg)
        end

        cfg.S.args = cfg.S.args[2:end]
        return QuantSymEx(cfg)

    elseif inst.head == :(=)

        if ~isa(inst.args[2], Expr)
            CAssign(cfg.σ, inst.args[1], inst.args[2])
            cfg.S.args = cfg.S.args[2:end]
            return QuantSymEx(cfg)
        end


        if inst.args[2].head == :comprehension # assume that comprehension doesn't contain conditionals
            inst.args[2].args[1].args[1] = postwalk(
                x -> x isa Symbol ? x == inst.args[2].args[1].args[2].args[1] ? Expr(:$, :($x)) : x : x,
                inst.args[2].args[1].args[1]
            )
            inst.args[2].args[1].args[1] = Expr(:block, Expr(:quote, inst.args[2].args[1].args[1]))
            qprogs = CEval(cfg.σ, inst.args[2])
            n_qprog = length(qprogs)
            temp = Vector{Union{Z3.ExprAllocated, Int}}(undef, n_qprog)
            S = copy(cfg.S)
            for j in 1:n_qprog
                cfg.S.args = [qprogs[j]]
                cfg = QuantSymEx(cfg)[1]
                temp[j] = cfg.σ[:__res__]
            end
            cfg.S.args = S.args[2:end]
            CAssign(cfg.σ, inst.args[1], temp)
            return QuantSymEx(cfg)
        elseif inst.args[2].head == :call # assume that subroutine doesn't contain conditionals
            S = copy(cfg.S)
            cfg.S.args = [inst.args[2]]
            cfg = QuantSymEx(cfg)[1]
            CAssign(cfg.σ, inst.args[1], cfg.σ[:__res__])
            cfg.S.args = S.args[2:end]
            return QuantSymEx(cfg)
        else
            CAssign(cfg.σ, inst.args[1], inst.args[2])
            cfg.S.args = cfg.S.args[2:end]
            return QuantSymEx(cfg)
        end

    elseif inst.head == :if
        ϕ = CEval(cfg.σ, inst.args[1])
        S1 = copy(cfg.S)
        S1.args = [inst.args[2].args;cfg.S.args[2:end]]
        S2 = copy(cfg.S)
        if length(inst.args) == 3
            S2.args = [inst.args[3].args;cfg.S.args[2:end]]
        else
            S2.args = cfg.S.args[2:end]
        end
        if ϕ isa Bool
            if ϕ
                cfg1 = SymConfig(S1, cfg.σ, cfg.ρ, cfg.P, cfg.ϕ, cfg.ctx)
                return QuantSymEx(cfg1)
            else
                cfg2 = SymConfig(S2, cfg.σ, cfg.ρ, cfg.P, cfg.ϕ, cfg.ctx)
                return QuantSymEx(cfg2)
            end
        end
        cfg1 = SymConfig(S1, cfg.σ, cfg.ρ, cfg.P, (cfg.ϕ[1] & ϕ, cfg.ϕ[2]), cfg.ctx)
        cfg2 = SymConfig(S2, cfg.σ, cfg.ρ, cfg.P, (cfg.ϕ[1] & not(ϕ), cfg.ϕ[2]), cfg.ctx)
        return vcat(QuantSymEx(cfg1)..., QuantSymEx(cfg2)...)
    elseif inst.head == :for
        inst.args[2] = postwalk(
            x -> x isa Symbol ? x == inst.args[1].args[1] ? Expr(:$, :($x)) : x : x,
            inst.args[2]
        )
        new_S = Expr(:block, Expr(:quote, inst.args[2]))
        qprogs = [CEval(CState(Dict([(inst.args[1].args[1], j)])), new_S) for j in CEval(cfg.σ, inst.args[1])]
        #qprogs = CEval(cfg.σ, Expr(:comprehension, new_S, inst.args[1]))
        n_qprog = length(qprogs)
        S = copy(cfg.S)
        for j in 1:n_qprog
            cfg.S.args = qprogs[j].args
            cfg = QuantSymEx(cfg)[1]
        end
        #S = copy(cfg.S)
        #cfg.S.args = [vcat([s.args for s in qprogs]...);cfg.S.args[2:end]]
        #@time QuantSymEx(cfg)
        cfg.S.args = S.args[2:end]
        return QuantSymEx(cfg)
    else
        cfg.σ[:__res__] = CEval(cfg.σ, inst)
        cfg.S.args = cfg.S.args[2:end]
        return QuantSymEx(cfg)
    end

end