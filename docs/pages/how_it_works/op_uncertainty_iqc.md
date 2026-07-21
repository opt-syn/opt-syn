---

tocdepth: 1

---


# Operators and Uncertainty




## Operator Classes

The individual operators $F_i$ in the inclusion problem may be contained in known operator classes. The optimization problem in {ref}:dr


The specific operator $F$ in the inclusion problem may not be a-priori known. The linear system $G$ should ensure that the algorithm $(F, G)$ satisfies the previous properties for all $F$ in a desired class of operators  $\F$. 

<!-- Example properties for functions $f_i$ include  convexity, strong convexity, and smoothness. Properties for operators $F_i$ include monotonicity, strong monotonicity,  -->

An instance of this operator class description are that $\F_1$ is the set of subdifferentials of  $m$-strongly convex and $L$-smooth function, and $F_2$ is the set of subdifferentials of a proper, closed, and convex functions. The class $\F$ is then set of operators $(F_1, F_2) \in (\F_1, \F_2)$ such that the zero-inclusion problem is feasible $(\exists \beta^* \mid 0 \in F_1(\beta^*) + F_2(\beta^*))$. 

A full listing of supported operator class descriptions is available in {doc}`Operators <../usage/problem_formulation/operators>`.




## Integral Quadratic Constraints

<!-- Membership of $F \in \F$ adds a restriction on the possible set of $(w, z)$ sequences obser  -->




## Computational Verification

 