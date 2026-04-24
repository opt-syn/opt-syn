# Operator 
Operators in the optimization algorithm.

Each operator F is a (possibly set-valued) map w \in F(z)

Sequences (w_k, z_k) satisfying w_k \in F(z_k) for all times k must obey a 
set of valid relations. An example of such a relation is the subgradient 
inequality for convex functions:

f(y) >= f(x) + <w', x-y> for all (x, y, w) with w \in \partial f(x)

These relations can be used to define valid Integral Quadratic Constraints 
that must be satisfied by sequences (w_k, z_k).

All operators support a coordinate-wise partition 'c': **TODO: finish this explanation for c, integrate into the feature**

## Supported Operators

The following operators can be described:

### General Set-Valued Maps

Operators w \in F(z) can satisfy combinations of the following constraints:
- Monotone
- mu-strongly monotone 
- mu-weakly monotone 
- beta-cocoercive
- Lipschitz
- Inverse Lipschitz

To describe an operator that is 1-strongly-monotone and 4-Lipschitz, use
`op = op_gen({'monotone', 1, 'lipschitz', 4})`.

### Subdifferentials

Operators w \in \partial f(z)  where there exist constants -inf < m <= L <= inf
such that the following functions are both proper, convex, and closed:

f(z) - m norm(z, 2)^2 and L norm(z, 2)^2 - f(z)

Convex functions have (m, L) = (0, inf). Instances of these convex functions
include indicator functions over convex sets.

The log-sum-exp function f(z) = log(sum(e^(b z)) has a subdifferential described by (0, b).


An operator characterized by (m, L) can be invoked by  `op = op_sml(m, L)`. To restrict the set of IQCs to causal (simpler computation, possible conservatism), use
`op = op_sml_causal(m, L)`.

### Quadratics

Quadratics are a specialization of the class (m, L).

The quadratic function f(z) = 1/2 (z-z0)' Q (z-z0) where Q is symmetric may be described by
m = min(eig(Q)) and L = max(eig(Q)).

Quadratics arise in least squares problems.

Prior knowledge of the quadratic nature of an operator can reduce conservatism as compared to using a general (m, L) description.

To invoke quadratics, use `op = op_quad(m, L)` or `op = op_quad_causal(m, L)`.

### Equality constraints

An equality constraint Ez = b can be imposed in the optimization problem.

It is assumed that the matrix E has full row rank

Equality constraints will be abstracted using the minimal and maximal singular values of E (smin, smax)

An equality constrained operator can be called by `op = op_eq(smin, smax)` or `op = op_eq_causal(smin, smax)`

Equality constraints can be used to analyze and generate primal-dual algorithms.
