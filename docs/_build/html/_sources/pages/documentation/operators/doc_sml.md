# Subdifferentials



These classes describe valid relations satisfied by  subdifferentials of functions in $S_{m, L}$.


The non-causal implementations are less conservative, but are more computationally intensive as compared to the causal implementation. 

The `order` supplied to Analysis is a single integer for causal (number of lags), and a pair of integers (number of primal lags, number of dual lags) for non-causal. Causal is equivalent a noncausal implementation with `order` = (number of lags, 0). 

## Non-causal

```{eval-rst}
.. mat:autoclass :: operator.op_sml  
    :members:
```

The set of p.c.c. functions is an instance of the class $\mathcal{S}_{0, \infty}$.
```{eval-rst}
.. mat:autoclass :: operator.op_pcc
    :members:
```

## Causal

```{eval-rst}
.. mat:autoclass :: operator.op_sml_causal  
    :members:
```

## Common Routines

```{eval-rst}
.. mat:autoclass :: operator.op_sml_interface  
    :members:
```