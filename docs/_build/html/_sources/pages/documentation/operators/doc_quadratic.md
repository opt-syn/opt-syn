# Gradients of Quadratics

These routines are appropriate for operators $F(z) = Q z - b$, where $Q$ is a symmetric matrix with eigenvalues between $m$ and $L$ and $b$ is a vector. $F$ is the gradient of a quadratic function $f(z) = \frac{1}{2}z^\top Q z - b^\top z + f_0$ for some constant $f_0$. 

A special case is the least-squares penalty $f(z) = \frac{1}{2} \norm{z-z_0}_2^2$, which is a quadratic with m=L.

## Noncausal


```{eval-rst}
.. mat:autoclass :: operator.op_quad  
    :members:
```

## Causal

```{eval-rst}
.. mat:autoclass :: operator.op_quad_causal  
    :members:
```


