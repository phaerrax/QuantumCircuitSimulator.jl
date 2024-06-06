# QuantumCircuitSimulator.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

With this Julia package you can load quantum circuits from OpenQASM files and
simulate them with tensor networks, using the
[ITensor](https://github.com/ITensor/ITensors.jl) library.

The package allows you to simulate noiseless or noisy circuits, and it also contains a rudimental implementation of the tensor-network error-mitigation technique developed in [1].


# References 
1. Sergei Filippov, Matea Leahy, Matteo A. C. Rossi, Guillermo García-Pérez. *Scalable tensor-network error mitigation for near-term quantum computing*. [arXiv:2307.11740](https://arxiv.org/abs/2307.11740).

