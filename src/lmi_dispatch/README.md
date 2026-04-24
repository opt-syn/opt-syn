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



