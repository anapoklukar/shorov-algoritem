# Shor's Algorithm

**Authors:** Blaž Grilj, Ana Poklukar

**Date:** November 2025

---

This project implements and analyzes **Shor's algorithm**, the quantum algorithm for factoring integers, using Qiskit. We built the full quantum circuit, including quantum phase estimation with a modular exponentiation oracle, ran it on both a local simulator and a real IBM quantum computer, and recovered the factors using continued fractions. We also compared the results against classical factoring methods to see how the quantum approach stacks up in practice.

The repository includes the full implementation along with our written report and presentation slides. Key components include:
- **Quantum Circuit:** QPE with a controlled modular-exponentiation oracle for period finding.
- **Execution:** local simulation via AerSimulator, or real hardware via IBM Quantum.
- **Classical Comparison:** brute-force search and Pollard's rho for reference.
- **Benchmarking:** success rate and runtime analysis over repeated trials.
