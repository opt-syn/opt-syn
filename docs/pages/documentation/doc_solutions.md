# Solution


The output of a {meth}`solve_single`, {meth}`bisect`, or {meth}`alternate` run from a {class}`opt_analysis` and {class}`opt_synthesis` {doc}`managers <doc_manager>` is a solution structure. The solution structure is of type {class}`opt_solution`, and is common among analysis and synthesis. Dedicated solution certificates for analysis and synthesis respectively are stored in the {attr}`cert` field of {class}`opt_solution`.

## Solution
```{eval-rst}
.. mat:autoclass :: manager.containers.opt_solution
    :members:
```

## Analysis certificate
```{eval-rst}
.. mat:autoclass :: manager.containers.cert_analysis
    :members:    
```
## Synthesis certificate

```{eval-rst}
.. mat:autoclass :: manager.containers.cert_synthesis 
    :members:
```