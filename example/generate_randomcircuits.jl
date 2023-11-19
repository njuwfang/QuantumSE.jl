ENV["JULIA_PKG_USING_AUTOINSTALL"] = "yes"
using StatsBase: sample

mkpath("./stim_benchmark/")

function generate_random_circuit(n)
    circuit = ""
    for j in 1:n
        circuit *= "H$(join(" $(k)" for k in sample(1:n,div(n,3)+1,replace=false)))\n"
        circuit *= "S$(join(" $(k)" for k in sample(1:n,div(n,3)+1,replace=false)))\n"
        circuit *= "CX$(join(" $(k)" for k in sample(1:n,10,replace=false)))\n"
        circuit *= "M$(join(" $(k)" for k in sample(1:n,div(n,20)+1,replace=false)))\n"
    end
    circuit *= "M$(join(" $(k)" for k in 1:n))"
    circuit
end

for j in 10:10:1000
    a = generate_random_circuit(j)
    open("./stim_benchmark/random$(j).stim", "w") do io
        println(io, a)
    end
end
