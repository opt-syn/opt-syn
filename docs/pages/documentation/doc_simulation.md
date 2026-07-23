# Simulation

 

Simulation routines are explored in the  {doc}`Simulation <../usage/simulation>` page.



## Simulation

`alg_sim` executes a given optimization algorithm, and stores the output in an `alg_sim_out` object.
```{eval-rst}
.. mat:autoclass :: simulator.alg_sim   
    :members:
```
```{eval-rst}
.. mat:autoclass :: simulator.alg_sim_out  
    :members:
```

<!--      -->
## Operators

The operators are defined by `op_sim` classes. Each operator $F$ has three core routines:


```{list-table}
:header-rows: 1

* - Evaluation
  - Name
  - Operation
* - Forward
  - `fw`
  - $z \mapsto F(z)$,
* - Backward
  - `bw`
  - $z \mapsto (I - D F)^{-1} (z)$
* - Function 
  - `f` 
  - $z \mapsto f(z)$
```
Function evaluation is supported if the operator $F$ is the subdifferential of a function $f$. When $F$ is the psuedogradient of a game with multiple agents, $f$ can be defined as the vector of payoff functions for each agent. If $f$ is undefined, then `f` returns the empty set `[]`. 



```{eval-rst}
.. mat:automodule :: simulator.op_sim
    :members:   
    :show-inheritance:
```