# Simulation

 

Simulation routines are explored in the  {doc}`Solve <../usage/simulation>` page.



## Simulate and Plots

`alg_sim` executes a given optimization algorithm. 
```{eval-rst}
.. mat:autoclass :: simulator.alg_sim   
    :members:
```


<!-- `alg_plotter` creates plots from an `alg_sim`-generated trajectory.

```{eval-rst}
.. mat:automodule :: simulator.alg_plotter        
``` -->


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