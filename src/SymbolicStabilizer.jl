using Z3
using LinearAlgebra: I

@inline _mod(j) = (j-1)>>6 + 1
@inline _rem(j) = UInt64(0x1) << ((j-1)&(0x3f))

_len2(j) = Int(ceil(log2(j)))

@inline _bv_val(ctx, v::Integer) = bv_val(ctx, v, 1)
@inline _bv_const(ctx, s::String) = bv_const(ctx, s, 1)

struct QState
    num_qubits::UInt32

    xzs::Matrix{UInt64}

    phases::Vector{Z3.ExprAllocated}

    ctx::Z3.ContextAllocated

    QState(num_qubits::Integer, Tableau::Matrix{Bool}, phases::Vector{Z3.ExprAllocated}, ctx::Z3.ContextAllocated) = begin
        @assert size(Tableau) == (2*num_qubits, 2*num_qubits)
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

    QState(num_qubits::Integer, Tableau::Matrix{Bool}, Phases::Vector{Bool}, ctx::Z3.ContextAllocated) = begin
        phases = Vector{Z3.ExprAllocated}(undef, 2*num_qubits)
        for j in 1:2*num_qubits
            phases[j] = _bv_val(ctx, 1*Phases[j])
        end

        QState(num_qubits, Tableau, phases, ctx)
    end

    QState(num_qubits::Integer, ctx::Z3.ContextAllocated) = begin
        Tableau = Matrix{Bool}(I, 2*num_qubits, 2*num_qubits)
        Phases = zeros(Bool, 2*num_qubits)
        
        QState(num_qubits, Tableau, Phases, ctx)
    end

    QState(q::QState) = begin
        new(q.num_qubits, copy(q.xzs), copy(q.phases), q.ctx)
    end

    QState(num_qubits::Integer, ctx::Z3.ContextAllocated) = begin
        Tableau = I + Zeros(Bool, 2*num_qubits, 2*num_qubits)
        phases = [_bv_val(ctx, 0) for j in 1:2*num_qubits]
        QState(num_qubits, Tableau, phases, ctx)
    end
end

Base.copy(q::QState) = QState(q)

function update!(q::QState, q0::QState)
    q.xzs[:] = q0.xzs
    q.phases[:] = q0.phases
end

@inline function rowset!(q::QState, i, b)
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

@inline function rowcopy!(q::QState, i, j)
    (i == j) && return
    q.phases[i] = q.phases[j]
    @inbounds @simd for k in 1:size(q.xzs,1)
        q.xzs[k,i] = q.xzs[k,j]
    end
    
    nothing
end

@inline function rowswap!(q::QState, i, j)
    (i == j) && return
    q.phases[i], q.phases[j] = q.phases[j], q.phases[i]
    @inbounds @simd for k in 1:size(q.xzs,1)
        q.xzs[k,i], q.xzs[k,j] = q.xzs[k,j], q.xzs[k,i]
    end
    
    nothing
end

function rowmult!(q::QState, i, j)
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

function _canonicalize_gott!(q::QState)
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
        
    result = QState(num_qubits, Tableau, phases, ctx)

    _canonicalize_gott!(result)

    result
end

function print_full_tableau(q::QState)
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

function CNOT!(q::QState, b, c)
    len = size(q.xzs, 1)÷2
    b6 = _mod(b)
    c6 = _mod(c)
    pwb = _rem(b)
    pwc = _rem(c)

    @inbounds @simd for j in 1:2*q.num_qubits
        x1, z1, x2, z2 = q.xzs[b6,j]&pwb, q.xzs[b6+len,j]&pwb, q.xzs[c6,j]&pwc, q.xzs[c6+len,j]&pwc
        if ~iszero(x1)
            q.xzs[c6,j] ⊻= pwc
        end
        if ~iszero(z2)
            q.xzs[b6+len,j] ⊻= pwb
        end
        if ~iszero( (x1 & z1 & x2 & z2)  | (x1 & z2 &~(z1|x2)) )
            q.phases[j] = q.phases[j] ⊻ _bv_val(q.ctx, 1)
        end
    end

    nothing
end

function H!(q::QState, b)
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

function S!(q::QState, b)
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

function X!(q::QState, b)
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

function Z!(q::QState, b)
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

function Y!(q::QState, b)
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

function sX!(q::QState, b, s)
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

function sZ!(q::QState, b, s)
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

function sY!(q::QState, b, s)
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

function M!(q::QState, b, sym_name::String)
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
        
        return res * 2
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

function inject_errors(q::QState, max_num_errors::Integer, error_type::String)
    terms = Vector{Z3.ExprAllocated}(undef, q.num_qubits)
    sGate = error_type == "X" ? sX : sZ
    for j in 1:q.num_qubits
        terms[j] = sGate(q, j, "$(error_type)_error_$(j)")
    end

    q.constraints = q.constraints & (sum( (x -> concat(bv_val(q.ctx, 0, _len2(q.num_qubits)), x)).(terms) ) <= bv_val(q.ctx, max_num_errors, _len2(q.num_qubits)+1))
    
    nothing
end

function inject_errors(q::QState, max_num_errors::Integer, error_type::String)
    terms = Vector{Z3.ExprAllocated}(undef, q.num_qubits)
    sGate = error_type == "X" ? sX : sZ
    e = [_bv_const(q.ctx, "$(error_type)_error_$(j)") for j in 1:q.num_qubits]
    for j in 1:q.num_qubits
        sGate(q, j, e[j])
    end

    return (sum( (x -> concat(bv_val(q.ctx, 0, _len2(q.num_qubits)), x)).(e) ) <= bv_val(q.ctx, max_num_errors, _len2(q.num_qubits)+1))
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

function check_state_equivalence(q1::QState, q2::QState, assumptions::Tuple{Z3.ExprAllocated, Z3.ExprAllocated, Z3.ExprAllocated}, slv_backend_cmd::Cmd=`yices-smt2`)
    q1.num_qubits == q2.num_qubits || error("The number of qubits does not match")
    
    ranges = q1.num_qubits+1:2*q1.num_qubits

    _equal(q1.xzs, q2.xzs, ranges) || error("The Stabilizer does not match, Your code is wrong even without error insertion")

    _canonicalize_gott!(q1)
    _canonicalize_gott!(q2)

    slv = Solver(q1.ctx)

    conjecture = reduce(&, [simplify(q1.phases[j] ⊻ q2.phases[j]) == _bv_val(q1.ctx, 0) for j in ranges])

    conjecture = assumptions[1] & (not(assumptions[2]) | (assumptions[3] & not(conjecture)))

    add(slv, conjecture)

    smt2_file_name = "_temp_check_equivalence_"

    open(smt2_file_name*".smt2", "w") do io
		println(io, "(set-logic QF_BV)")
		println(io, slv)
		println(io, "(check-sat)")
		println(io, "(get-model)")
		println(io, "(exit)")
	end

    run(pipeline(
        slv_backend_cmd,
        stdin=smt2_file_name*".smt2",
        stdout=smt2_file_name*".output",
        append=true
    ))

    nothing
end