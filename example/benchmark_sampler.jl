using Pkg
Pkg.activate("../")

using QuantumSE
using Z3

ctx = Z3.Context()

parse_stim("/home/fangw/repos/stim_benchmark/random10.stim", ctx)

open("benchmark_quantumse.dat", "w") do io
  println(io, "nq init_time sample_time")
  for j in 10:10:1000
    t0 = time()
    rs = parse_stim("/home/fangw/repos/stim_benchmark/random$(j).stim", ctx)
    sampler = Sampler(rs, 0.5)
    t1 = time()
    sampler(1000)
    t2 = time()
    println(io, "$(j) $(t1-t0) $(t2-t1)")
    println("$(j) $(t1-t0) $(t2-t1)")
  end
end
