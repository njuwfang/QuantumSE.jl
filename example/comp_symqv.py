import time

import numpy as np

from symqv.lib.expressions.qbit import Qbits
from symqv.lib.models.circuit import Circuit, Method
from symqv.lib.operations.gates import H, Z, X
from symqv.lib.solver import SpecificationType

def prove_id(n: int, delta=0.0001):
    # Initialize circuit
    qbits = Qbits([f'q{i}' for i in range(n)])

    circuit = Circuit(qbits,
                      [
                          [H(qbit) for qbit in qbits],
                          [X(qbit) for qbit in qbits],
                          [X(qbit) for qbit in qbits],
                          [H(qbit) for qbit in qbits],
                      ],
                      delta=delta)

    # Symbolic execution
    final_qbits = circuit.get_final_qbits()

    circuit.set_specification(
        [(final_qbits[i], qbits[i]) for i
         in
         range(n)],
        SpecificationType.equality_pair_list)

    sat,_,_ = circuit.prove(method=Method.qbit_sequence_model,
                  overapproximation=True)

    print(sat)

if __name__ == "__main__":
    for i in range(1,21):
        times = []

        for _ in range(1):
            start = time.time()
            prove_id(i)
            times.append(time.time() - start)

        print(f'Runtime for {i}:', np.mean(times))
