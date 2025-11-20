#!/usr/bin/env python3

import argparse
import random
import numpy as np
from math import gcd, floor
from fractions import Fraction
from IPython.display import display

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

    return qc.to_gate(label=f"{a}^(2^{binary_power}) mod {N}").control()


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
    #print(f"Measured {measurement} → {decimal_value} → r ≈ {r}")
    if not isinstance(r, int) or r % 2 != 0:
        return None

    x = pow(a, r // 2, N)
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
def factor_with_shor(N, shots=1024, use_ibmq=False):
    n = floor(np.log2(N - 1)) + 1

    if use_ibmq:
        print("\nUsing IBM backend...")
        service = QiskitRuntimeService()
        backend = service.backend("ibm_fez")

        sampler = Sampler(backend)
        sampler.options.dynamical_decoupling.enable = True
        sampler.options.dynamical_decoupling.sequence_type = "XpXm"
        sampler.options.twirling.enable_gates = True

    # Choose random a coprime to N
    while True:
        a = random.randint(2, N - 1)
        if gcd(a, N) == 1:
            break

    print(f"\nChosen a = {a}")

    # Build QPE circuit
    qpe = build_qpe_circuit(a, N)
    try:
        qpe.draw(output="mpl", fold=-1).savefig("shor_qpe_circuit.png")
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

    plot_histogram(counts, figsize=(35, 5)).savefig("shor_measurement_histogram.png")
    sorted_measurements = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    for meas, _ in sorted_measurements:
        factors = try_extract_factors(a, meas, n, N)
        if factors:
            print(f"\nSUCCESS: {N} = {factors[0]} x {factors[1]}")
            return factors

    print("\nFailed to factor N.")
    return None


# ----------------------------------------------------------
# Command-line argument handling
# ----------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Run Shor's Algorithm.")
    parser.add_argument("-N", type=int, required=True, help="Number to factor")
    parser.add_argument("-s", "--shots", type=int, default=1024, help="Number of shots")
    parser.add_argument("--ibm", action="store_true", help="Use IBM backend instead of local simulator")

    args = parser.parse_args()

    factor_with_shor(N=args.N, shots=args.shots, use_ibmq=args.ibm)


if __name__ == "__main__":
    main()
