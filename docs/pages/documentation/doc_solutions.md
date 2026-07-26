# Solutions


The output of a {meth}`solve_single`, {meth}`bisect`, or {meth}`alternate` run from a {class}`opt_analysis` and {class}`opt_synthesis` {doc}`managers <doc_manager>` is a solution structure. The solution structure is of type {class}`opt_solution`, and is common among Analysis and Synthesis. Dedicated solution certificates for Analysis and Synthesis respectively are stored in the {attr}`cert` field of {class}`opt_solution`.

## Solution
```{eval-rst}
.. mat:autoclass :: manager.containers.opt_solution
    :members:
```

## Analysis Certificate
```{eval-rst}
.. mat:autoclass :: manager.containers.cert_analysis
    :members:    
```
## Synthesis Certificate

```{eval-rst}
.. mat:autoclass :: manager.containers.cert_synthesis 
    :members:
```