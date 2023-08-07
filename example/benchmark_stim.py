import stim
from time import time

circuit = stim.Circuit.from_file(f"./stim_benchmark/random10.stim")
sampler = circuit.compile_sampler()
sampler.sample(shots=10000)

with open("./benchmark_stim.dat", "w") as io:
    print("nq init_time sample_time", file=io)
    for j in range(10,1001,10):
        t0 = time()
        circuit = stim.Circuit.from_file(f"./stim_benchmark/random{j}.stim")
        sampler = circuit.compile_sampler()
        t1 = time()
        sampler.sample(shots=10000)
        t2 = time()
        print(f"{j} {t1-t0} {t2-t1}", file=io)
        print(f"{j} {t1-t0} {t2-t1}")

