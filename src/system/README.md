# SRC

All codes for the opt-syn project: analysis and synthesis of optimization algorithms.

An optimization algorithm is a procedure that finds a point beta solving the fixed-point equation

0 \in \sum_{i=1}^s F_i(\beta)


These fixed point equations arise in applications including optimization, control, operations research, imaging.

The maps F could arise from subgradients of convex functions, monotone operators, or equality constraints

## Setting up the problem

An algorithmic interconnection is described by five attributes:
1. Operators
2. Network
3. Controller
4. Performance (optional)
5. Tracking (optional)

The operators are the oracles F in the optimization problem.


## Operator

## LMI_DISPATCH: 

contains the lmi routines for analysis and synthesis


# System types (from opt_system):

Implemented:
- Linear Time Invariant (LTI)
- Switched (robust)
- Switched (stochastic: Markov Jump Linear System)

Periodic is a special case of Switched Robust

Future:
- LPV (polytopic)
- LPV (linear fractional)

# Performance Specifications (from spec):

Implemented:
- stability:      exponential stability of interconnection
- quad:           quadratic performance
- e2e:            energy to energy (l2) gain

Future:
- h2:             primal H2 (maybe dual as well?)
- p2p:            peak to peak gain (l2)
- e2p:            energy to peak  gain (generalized h2)
- e2p:            energy to peak  gain (generalized h2)


# Implemented Pairs

| **System Type**         | **Analysis** | **Synthesis** |
|-------------------------|--------------|---------------|
| LTI                     |              |               |
| Stochastic (robust)     |              |               |
| Switched (stochastic)   |              |               |
| LPV (polytopic)         |              |               |
| LPV (linear fractional) |              |               |
