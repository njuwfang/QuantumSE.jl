module QSE

include("QuantumProgram.jl")
export @qprog, QProg, QEmpty

include("CSSCode.jl")
export GF2

include("SymbolicStabilizer.jl")
export M, X, Y, Z, sX, sY, sZ, S, H, CNOT,
    QState, from_stabilizer, print_full_tableau, update!, inject_errors, from_css_code,
    _bv_val, _bv_const, _len2, check_state_equivalence, _sum

include("SymbolicExecution.jl")
export CState, SymConfig, QuantSymEx, CEval

include("LinearGroup.jl")
export AbstractGroup, AbstraLinearGroup, L2, PL2, pgl2, psl2, jacobi4squares

end # module QSE
