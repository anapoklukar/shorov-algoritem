#!/usr/bin/env python3

import argparse
import random
import numpy as np
from math import gcd, floor
from fractions import Fraction
from IPython.display import display
import time   ### ADDED

from qiskit_aer import AerSimulator
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, Sampler
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import QuantumRegister, ClassicalRegister
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram

## Uncomment and set up your IBM Quantum account token if you plan to use IBM backends
# QiskitRuntimeService.save_account(
#     token="API-TOKEN-HERE",
#     set_as_default = True,
#     overwrite= True
# )


# ----------------------------------------------------------
# Helper: Continued fraction denominator
# ----------------------------------------------------------
def denominator(decimal_value, n, N):
    return (Fraction(decimal_value / 2 ** (2 * n)).limit_denominator(N)).denominator


# ----------------------------------------------------------
# Controlled multiply gate a^(2^k) mod N
# ----------------------------------------------------------
def ctrl_mult_gate(a, binary_power, N):
    n = N.bit_length()
    value = pow(a, 2 ** binary_power, N)
    qc = QuantumCircuit(n)

    for i in range(n):
        if (value >> i) & 1:
            qc.x(i)

    return qc.to_gate(label=f"{a}^{2**binary_power} mod {N}").control()


# ----------------------------------------------------------
# Apply U_f controlled-multiplication
# ----------------------------------------------------------
def apply_Uf(circ, qreg, treg, a, N):
    n = N.bit_length()
    for k in range(2 * n):
        gate = ctrl_mult_gate(a, k, N)
        circ.append(gate, [qreg[k]] + list(treg))


# ----------------------------------------------------------
# Build QPE circuit to estimate period r
# ----------------------------------------------------------
def build_qpe_circuit(a, N):
    n = N.bit_length()
    q = QuantumRegister(2 * n, "q")
    t = QuantumRegister(n, "t")
    c = ClassicalRegister(2 * n, "c")

    qc = QuantumCircuit(q, t, c)

    qc.x(t[0])      # initialize |1>
    qc.h(q)         # Hadamards on counting register
    qc.barrier()

    apply_Uf(qc, q, t, a, N)
    qc.barrier()

    qc.append(QFT(2 * n, inverse=True), q)
    qc.measure(q, c)
    return qc


# ----------------------------------------------------------
# Extract factors from measurement
# ----------------------------------------------------------
def try_extract_factors(a, measurement, n, N):
    decimal_value = int(measurement, 2)
    r = denominator(decimal_value, n, N)
    #print(f"Trying measurement {measurement} (decimal {decimal_value}) → r = {r}")
    #print(f"Measured {measurement} → {decimal_value} → r ≈ {r}")
    if not isinstance(r, int) or r % 2 != 0:
        return None

    x = pow(a, r // 2, N)
    #print(f"Computed x = a^(r/2) mod N = {x}")
    f1 = gcd(x - 1, N)
    f2 = gcd(x + 1, N)

    if f1 not in (1, N) and f2 not in (1, N) and f1 * f2 == N:
        print("Period r found:", r)
        return (f1, f2)

    if f1 not in (1, N) and N % f1 == 0:
        print("Period r found:", r)
        return (f1, N // f1)

    if f2 not in (1, N) and N % f2 == 0:
        print("Period r found:", r)
        return (f2, N // f2)

    return None


# ----------------------------------------------------------
# Shor algorithm main function
# ----------------------------------------------------------
def factor_with_shor(N, shots=1024, use_ibmq=False, draw=False):   ### ADDED draw argument
    n = floor(np.log2(N - 1)) + 1

    if use_ibmq:
        print("\nUsing IBM backend...")
        service = QiskitRuntimeService()
        backend = service.backend("ibm_torino")

        sampler = Sampler(backend)
        sampler.options.dynamical_decoupling.enable = True
        sampler.options.dynamical_decoupling.sequence_type = "XpXm"
        sampler.options.twirling.enable_gates = True

    # Choose random a coprime to N
    a = random.randint(2, N - 1)

    while True:
        d = gcd(a, N)
        if d != 1:
            # nontrivial factor found
            print(f"Lucky guess of a = {a}, found factor {d}")
            print(f"\nSUCCESS: {N} = {d} x {N // d}")
            return (d, N // d)

        # gcd == 1 → good a
        break
    
    # Build QPE circuit
    qpe = build_qpe_circuit(a, N)
    #print("Number of qubits used:", qpe.num_qubits)
    ### ADDED drawing
    if draw:
        try:
            qpe.draw(output="mpl", fold=-1).savefig("shor_qpe_circuit.pdf", format="pdf")
        except:
            pass

    # Run
    if use_ibmq:
        pm = generate_preset_pass_manager(optimization_level=2, backend=backend)
        transpiled = pm.run(qpe)
        job = sampler.run([transpiled], shots=shots)
        counts = job.result()[0].data["c"].get_counts()
    else:
        sim = AerSimulator()
        compiled = transpile(qpe, sim)
        result = sim.run(compiled, shots=shots).result()
        counts = result.get_counts()

    if draw:
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(12, 5))         # smaller & tighter
            plot_histogram(counts, ax=ax)

            ax.set_ylabel("Število meritev")

            fig.tight_layout(pad=0.4)                       # tighter layout
            fig.savefig("shor_histogram.pdf", format="pdf")
        except:
            pass

    # Dictionary of bitstrings and their counts to keep
    counts_keep = {}
    # Threshold to filter
    threshold = np.max(list(counts.values())) / 2
    
    for key, value in counts.items():
        if value > threshold:
            counts_keep[key] = value
    
    #print(counts_keep)

    factors = try_extract_factors(a, list(counts_keep.keys())[0], n, N)

    if factors:
        print(f"\nSUCCESS: {N} = {factors[0]} x {factors[1]}")
        return factors

    print("\nFailed to factor N from a single measurement.")
    return None


# ----------------------------------------------------------
# BENCHMARK (ADDED)
# ----------------------------------------------------------
def benchmark(N, shots=1024, runs=100, use_ibmq=False):
    print(f"\nRunning benchmark: {runs} runs for N={N}")
    successes = 0
    failures = 0
    total_time = 0.0

    for _ in range(runs):
        start = time.time()
        result = factor_with_shor(N, shots=shots, use_ibmq=use_ibmq, draw=False)
        total_time += time.time() - start

        if result == (3, 5) or result == (5, 3):
            successes += 1
        else:
            failures += 1

    avg_time = total_time / runs

    print("\n--- Benchmark Results ---")
    print(f"Successes: {successes}")
    print(f"Failures: {failures}")
    print(f"Success rate: {successes/runs*100:.2f}%")
    print(f"Average time: {avg_time:.4f} seconds")
    print("-------------------------\n")


# ----------------------------------------------------------
# Command-line argument handling
# ----------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Run Shor's Algorithm.")
    parser.add_argument("-N", type=int, required=True, help="Number to factor")
    parser.add_argument("-s", "--shots", type=int, default=1024, help="Number of shots")
    parser.add_argument("--ibm", action="store_true", help="Use IBM backend instead of local simulator")

    parser.add_argument("--draw", action="store_true", help="Draw circuit & histogram")   ### ADDED
    parser.add_argument("--benchmark", type=int, nargs="?", const=100,
                        help="Benchmark mode (default 100 runs)")   ### ADDED

    args = parser.parse_args()

    ### ADDED logic
    if args.benchmark:
        benchmark(args.N, shots=args.shots, runs=args.benchmark, use_ibmq=args.ibm)
        return

    factor_with_shor(args.N, shots=args.shots, use_ibmq=args.ibm, draw=args.draw)   


if __name__ == "__main__":
    main()
