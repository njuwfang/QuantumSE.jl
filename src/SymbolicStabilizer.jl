using Z3
using LinearAlgebra: I, nullspace, rank
using SimpleGF2: rref
using InvertedIndices

@inline _mod(j) = (j-1)>>6 + 1
@inline _rem(j) = UInt64(0x1) << ((j-1)&(0x3f))

_len2(j) = Int(ceil(log2(j))) + 1

@inline _bv_val(ctx, v::Integer) = bv_val(ctx, v, 1)
@inline _bv_const(ctx, s::String) = bv_const(ctx, s, 1)

_sum(ctx, e, num_qubits) = sum( (x -> concat(bv_val(ctx, 0, _len2(num_qubits)), x)).(e) )

struct SymStabilizerState <: AbstractSymQuantumState
    num_qubits::UInt32

    xzs::Matrix{UInt64}

    phases::Vector{Z3.ExprAllocated}

    ctx::Z3.ContextAllocated

    SymStabilizerState(num_qubits::Integer, Tableau::Matrix{Bool}, phases::Vector{Z3.ExprAllocated}, ctx::Z3.ContextAllocated) = begin
        @assert size(Tableau) == (2*num_qubits, 2*num_qubits)
        println("Expected Len(phases) = $(2*num_qubits),; Got $(length(phases))")
        @assert length(phases) == 2*num_qubits
        
        len = num_qubits>>6 + 1
        xzs = zeros(UInt64, 2*len, 2*num_qubits+1)
        
        @simd for i in 1:2*num_qubits
            for j in 1:num_qubits
                j6 = _mod(j)
                pw = _rem(j)
                if ~iszero(Tableau[i,j])
                    xzs[j6,i] |= pw
                end
                if ~iszero(Tableau[i,j+num_qubits])
                    xzs[j6+len,i] |= pw
                end
            end
        end

        new(num_qubits, xzs, vcat(phases, _bv_val(ctx, 0)), ctx)
    end

    SymStabilizerState(num_qubits::Integer, Tableau::Matrix{Bool}, Phases::Vector{Bool}, ctx::Z3.ContextAllocated) = begin
        phases = Vector{Z3.ExprAllocated}(undef, 2*num_qubits)
        for j in 1:2*num_qubits
            phases[j] = _bv_val(ctx, 1*Phases[j])
        end

        SymStabilizerState(num_qubits, Tableau, phases, ctx)
    end

    #=SymStabilizerState(num_qubits::Integer, ctx::Z3.ContextAllocated) = begin
        Tableau = Matrix{Bool}(I, 2*num_qubits, 2*num_qubits)
        Phases = zeros(Bool, 2*num_qubits)
        
        SymStabilizerState(num_qubits, Tableau, Phases, ctx)
    end=#

    SymStabilizerState(q::SymStabilizerState) = begin
        new(q.num_qubits, copy(q.xzs), copy(q.phases), q.ctx)
    end

    SymStabilizerState(num_qubits::Integer, ctx::Z3.ContextAllocated) = begin
        Tableau = Matrix{Bool}(I, 2*num_qubits, 2*num_qubits)
        phases = [_bv_val(ctx, 0) for j in 1:2*num_qubits]
        SymStabilizerState(num_qubits, Tableau, phases, ctx)
    end
end

Base.copy(q::SymStabilizerState) = SymStabilizerState(q)

function update!(q::SymStabilizerState, q0::SymStabilizerState)
    q.xzs[:] = q0.xzs
    q.phases[:] = q0.phases
end

@inline function rowset!(q::SymStabilizerState, i, b)
    len = size(q.xzs, 1)÷2

    @inbounds @simd for j in 1:size(q.xzs, 1)
        q.xzs[j,i] = 0
    end
    q.phases[i] = _bv_val(q.ctx, 0)

    if b <= q.num_qubits
        q.xzs[_mod(b),i] = _rem(b)
    else
        q.xzs[len+_mod(b-q.num_qubits),i] = _rem(b-q.num_qubits)
    end

    nothing
end

@inline function rowcopy!(q::SymStabilizerState, i, j)
    (i == j) && return
    q.phases[i] = q.phases[j]
    @inbounds @simd for k in 1:size(q.xzs,1)
        q.xzs[k,i] = q.xzs[k,j]
    end
    
    nothing
end

@inline function rowswap!(q::SymStabilizerState, i, j)
    (i == j) && return
    q.phases[i], q.phases[j] = q.phases[j], q.phases[i]
    @inbounds @simd for k in 1:size(q.xzs,1)
        q.xzs[k,i], q.xzs[k,j] = q.xzs[k,j], q.xzs[k,i]
    end
    
    nothing
end

function rowmult!(q::SymStabilizerState, i, j)
    len = size(q.xzs, 1)÷2
    r = q.xzs[:,i]
    l = q.xzs[:,j]

    cnt1 = zero(q.xzs[1,1])
    cnt2 = zero(q.xzs[1,1])

    @inbounds @simd for k in 1:len
        x1, x2, z1, z2 = l[k], r[k], l[k+len], r[k+len]
        q.xzs[k,i] = newx1 = x1 ⊻ x2
        q.xzs[k+len,i] = newz1 = z1 ⊻ z2
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2 ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        cnt1 ⊻= anti_comm
    end

    extra_phase = ((count_ones(cnt1)⊻(count_ones(cnt2)<<1))&0x3)

    q.phases[i] = simplify(_bv_val(q.ctx, extra_phase÷2) ⊻ q.phases[i] ⊻ q.phases[j])

    nothing
end

function _canonicalize_gott!(q::SymStabilizerState)
    len = size(q.xzs, 1)÷2
    rest_j = UInt32[]
    perms = Tuple{UInt32,UInt32}[]

    @inbounds @simd for i in 1:q.num_qubits
        for j in 1:len
            q.xzs[j,i] = 0
            q.xzs[j+len,i] = 0
        end
        q.phases[i] = _bv_val(q.ctx, 0)
    end

    i = q.num_qubits + 1
    for j in 1:q.num_qubits
        j6 = _mod(j)
        pwj = _rem(j)
        k = findfirst(ii-> q.xzs[j6,ii]&pwj != 0, i:2*q.num_qubits)
        if k !== nothing
            k = k + i - 1
            rowswap!(q, k, i)
            push!(perms, (k, i))
            for k2 in vcat(q.num_qubits+1:i-1, i+1:2*q.num_qubits)
                if ~iszero(q.xzs[j6,k2]&pwj)
                    rowmult!(q, k2, i)
                end
            end
            q.xzs[j6+len,i-q.num_qubits] ⊻= pwj
            i += 1
        else
            push!(rest_j, j)
        end
    end

    for j in rest_j
        j6 = _mod(j)
        pwj = _rem(j)
        k = findfirst(ii-> q.xzs[j6+len,ii]&pwj != 0, i:2*q.num_qubits)
        if k !== nothing
            k = k + i - 1
            rowswap!(q, k, i)
            push!(perms, (k, i))
            for k2 in vcat(q.num_qubits+1:i-1, i+1:2*q.num_qubits)
                if ~iszero(q.xzs[j6+len,k2]&pwj)
                    rowmult!(q, k2, i)
                end
            end
            q.xzs[j6,i-q.num_qubits] ⊻= pwj
            i += 1
        end
    end

    for (k, i) in reverse(perms)
        rowswap!(q, k, i)
        rowswap!(q, k-q.num_qubits, i-q.num_qubits)
    end

    nothing
end

function from_stabilizer(num_qubits::Integer, Stabilizer::Matrix{Bool}, phases1::Vector{Z3.ExprAllocated}, ctx::Z3.ContextAllocated)
    Tableau = vcat(zeros(Bool, num_qubits, 2*num_qubits), Stabilizer)
    phases0 = Vector{Z3.ExprAllocated}(undef, num_qubits)
    for j in 1:num_qubits
        phases0[j] = _bv_val(ctx, 0)
    end
    phases = vcat(phases0, phases1)
        
    result = SymStabilizerState(num_qubits, Tableau, phases, ctx)

    _canonicalize_gott!(result)

    result
end

function logical_operators(H1, H2)
    n = size(H1,2)
    X = rref(H1)
    nx = rank(X)
    X_idxs = [findfirst(!iszero, X[j,:]) for j in 1:nx]
    Z = rref(H2[:,Not(X_idxs)])
    nz = rank(Z)
    Z_idxs = [[1:n...][Not(X_idxs)][findfirst(!iszero, Z[j,:])] for j in 1:nz]
    nl = n - nx - nz
    L = zeros(GF2, nl, n)
    L[1:nl, X_idxs] = X[1:nx,Not([X_idxs;Z_idxs])]'
    L[1:nl, Not([X_idxs;Z_idxs])] += I
    
    L
end

function from_css_code(HX, HZ, ctx::Z3.ContextAllocated)
    n = size(HX, 2)

    X = rref(HX)
    nx = rank(X)
    X = X[1:nx,:]

    Z = rref(HZ)
    nz = rank(Z)
    Z = Z[1:nz,:]

    LZ = logical_operators(X, Z)
    LX = logical_operators(Z, X)
    nl = size(LZ, 1)
    dx = minimum([length(findall(!iszero, LX[j,:])) for j in 1:nl])
    dz = minimum([length(findall(!iszero, LZ[j,:])) for j in 1:nl])
    println("[[n, k, dx, dz]] = [[$(n), $(nl), dx<$(dx), dz<$(dz)]]")
    
    stabilizer1 = Matrix{Bool}([[X;LX] zeros(GF2, nx+nl, n);zeros(GF2, nz, n) Z])
    phases1 = [_bv_val(ctx, 0) for j in 1:n]
    for j in 1:nl
        phases1[nx+j] = _bv_const(ctx, "lx$(j)")
    end
    ρ1 = from_stabilizer(n, stabilizer1, phases1, ctx)

    stabilizer2 = Matrix{Bool}([X zeros(GF2, nx, n);zeros(GF2, nz+nl, n) [Z;LZ]])
    phases2 = [_bv_val(ctx, 0) for j in 1:n]
    for j in 1:nl
        phases2[nx+nz+j] = _bv_const(ctx, "lz$(j)")
    end
    ρ2 = from_stabilizer(n, stabilizer2, phases2, ctx)

    ρ1, ρ2, dx, dz
end

function print_full_tableau(q::SymStabilizerState)
    len = size(q.xzs, 1)÷2
    
    for i in 1:2*q.num_qubits
        for j in 1:q.num_qubits
            j6 = _mod(j)
            pwj = _rem(j)
            print((q.xzs[j6,i]&pwj)>>((j-1)&(0x3f)))
        end
        for j in 1:q.num_qubits
            j6 = _mod(j)
            pwj = _rem(j)
            print((q.xzs[j6+len,i]&pwj)>>((j-1)&(0x3f)))
        end
        println(" | ", simplify(q.phases[i]))
        if i == q.num_qubits
            print("\n")
        end
    end
    nothing
end

function CNOT!(q::SymStabilizerState, b, c)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    c6 = _mod(c)
    pwb = _rem(b)
    pwc = _rem(c)

    @inbounds @simd for j in 1:2*q.num_qubits
        x1, z1, x2, z2 = q.xzs[b6,j]&pwb!=0, q.xzs[b6+len,j]&pwb!=0, q.xzs[c6,j]&pwc!=0, q.xzs[c6+len,j]&pwc!=0
        if x1
            q.xzs[c6,j] ⊻= pwc
        end
        if z2
            q.xzs[b6+len,j] ⊻= pwb
        end
        if (x1 && z1 && x2 && z2)  || (x1 && z2 && ~(z1||x2))
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end

    nothing
end

function H!(q::SymStabilizerState, b)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)
    
    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        q.xzs[b6,j] ⊻= x ⊻ z
        q.xzs[b6+len,j] ⊻= x ⊻ z
        if x!=0 && z!=0
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end

    nothing
end

function S!(q::SymStabilizerState, b)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)

    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        q.xzs[b6+len,j] ⊻= x
        if x!=0 && z!=0
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end
    
    nothing
end

function X!(q::SymStabilizerState, b)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)

    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        if ~iszero(z)
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end
    
    nothing
end

function Z!(q::SymStabilizerState, b)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)

    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        if ~iszero(x)
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end
    
    nothing
end

function Y!(q::SymStabilizerState, b)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)

    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        if ~iszero(x⊻z)
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end

    nothing
end

function sX!(q::SymStabilizerState, b, s)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)
    
    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        if ~iszero(z)
            q.phases[j] = q.phases[j] ⊻ s
        end
    end
    
    nothing
end

function sZ!(q::SymStabilizerState, b, s)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)
    
    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        if ~iszero(x)
            q.phases[j] = q.phases[j] ⊻ s
        end
    end
    
    nothing
end

function sY!(q::SymStabilizerState, b, s)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    pw = _rem(b)

    @inbounds @simd for j in 1:2*q.num_qubits
        x, z = q.xzs[b6,j]&pw, q.xzs[b6+len,j]&pw
        if ~iszero(x⊻z)
            q.phases[j] = q.phases[j] ⊻ s
        end
    end
    
    nothing
end

function M!(q::SymStabilizerState, b, sym_name::String)
    b6 = _mod(b)
    pw = _rem(b)
    p = findfirst(ii-> q.xzs[b6,ii+q.num_qubits]&pw != 0, 1:q.num_qubits)
    if p !== nothing
        rowcopy!(q, p, p+q.num_qubits)
        rowset!(q, p+q.num_qubits, b+q.num_qubits)
        
        res = _bv_const(q.ctx, sym_name)
        
        q.phases[p+q.num_qubits] = res

        @inbounds @simd for j in vcat(1:p-1, p+1:2*q.num_qubits)
            if ~iszero(q.xzs[b6,j]&pw)
                rowmult!(q, j, p)
            end
        end
        
        return res
    else
        m = findfirst(ii-> q.xzs[b6,ii]&pw != 0, 1:q.num_qubits)
        rowcopy!(q, 2*q.num_qubits+1, m+q.num_qubits)

        @inbounds for j in m+1:q.num_qubits
            if ~iszero(q.xzs[b6,j]&pw)
                rowmult!(q, 2*q.num_qubits+1, j+q.num_qubits)
            end
        end

        return q.phases[2*q.num_qubits+1]
    end
end

CNOT, H, S, X, Y, Z, sX, sY, sZ, M = CNOT!, H!, S!, X!, Y!, Z!, sX!, sY!, sZ!, M!

function inject_errors(q::SymStabilizerState, max_num_errors::Integer, error_type::String)
    terms = Vector{Z3.ExprAllocated}(undef, q.num_qubits)
    sGate = error_type == "X" ? sX : sZ
    for j in 1:q.num_qubits
        terms[j] = sGate(q, j, "$(error_type)_error_$(j)")
    end

    q.constraints = q.constraints & (sum( (x -> concat(bv_val(q.ctx, 0, _len2(q.num_qubits)), x)).(terms) ) <= bv_val(q.ctx, max_num_errors, _len2(q.num_qubits)+1))
    
    nothing
end

function inject_errors(q::SymStabilizerState, error_type::String)
    terms = Vector{Z3.ExprAllocated}(undef, q.num_qubits)
    sGate = error_type == "X" ? sX : sZ
    errors = [_bv_const(q.ctx, "$(error_type)_error_$(j)") for j in 1:q.num_qubits]
    for j in 1:q.num_qubits
        sGate(q, j, errors[j])
    end

    return errors
end

@inline function _equal(a, b, ranges)
    e = true
    for j in ranges
        for i in axes(a, 1)
            e &= a[i,j] == b[i,j]
        end
    end
    e
end

function check_state_equivalence(q1::SymStabilizerState, q2::SymStabilizerState, assumptions::Tuple{Z3.ExprAllocated, Z3.ExprAllocated, Z3.ExprAllocated}, slv_backend_cmd::Cmd=`bitwuzla -e 0 -SE kissat`;use_z3=false)
    if q1.num_qubits != q2.num_qubits
        @info "The number of qubits does not match"
        return false
    end
    
    ranges = q1.num_qubits+1:2*q1.num_qubits

    _canonicalize_gott!(q1)
    _canonicalize_gott!(q2)

    if ~_equal(q1.xzs, q2.xzs, ranges)
        @error "The Stabilizer does not match, the program is wrong even without error insertion"
        return false
    end

    slv = Solver(q1.ctx)

    add(slv, assumptions[1])
    add(slv, not(assumptions[2]))

    if use_z3
        res = check(slv)
        if res == Z3.sat
            return false
        end
    else

        smt2_file_name = "_temp_check_preconditons_"

        open(smt2_file_name*".smt2", "w") do io
	    	println(io, "(set-logic QF_BV)")
            println(io, "(set-option :produce-models true)")
	    	println(io, slv)
	    	println(io, "(check-sat)")
	    	println(io, "(get-model)")
	    	println(io, "(exit)")
	    end

        res_string = read(pipeline(`$(slv_backend_cmd) $(smt2_file_name*".smt2")`), String)

        if ~occursin("unsat", res_string)
            @info "The preconditions of external programs are not satisfied"
            open(smt2_file_name*".output", "w") do io
                println(io, res_string)
            end
            @error "The assignment that generates the bug has been written to ./$(smt2_file_name).output"
            return false
        end
    end

    Z3.reset(slv)

    println("Slv after pre-cond reset: $(slv)")

    conjecture = reduce(&, [simplify(q1.phases[j] ⊻ q2.phases[j]) == _bv_val(q1.ctx, 0) for j in ranges])
    
    println("Conjecture by reduce: $(conjecture)")
    conjecture = assumptions[1] & assumptions[3] & not(conjecture)

    println("A1: $(assumptions[1])")
    add(slv, conjecture)
    println("A3: $(assumptions[3])")

    if use_z3
        res = check(slv)
        if res == Z3.sat
            return false
        end
    else
        smt2_file_name = "_temp_check_equivalence_"

        open(smt2_file_name*".smt2", "w") do io
	    	println(io, "(set-logic QF_BV)")
            println(io, "(set-option :produce-models true)")
	    	println(io, slv)
	    	println(io, "(check-sat)")
	    	println(io, "(get-model)")
	    	println(io, "(exit)")
	    end

        res_string = read(pipeline(`$(slv_backend_cmd) $(smt2_file_name*".smt2")`), String)

        if ~occursin("unsat", res_string)
            @info "There exist some allowed errors that the program cannot correct"
            open(smt2_file_name*".output", "w") do io
                println(io, res_string)
            end
            @error "The assignment that generated the bug has been written to ./$(smt2_file_name).output"
            return false
        end
    end

    return true
end
