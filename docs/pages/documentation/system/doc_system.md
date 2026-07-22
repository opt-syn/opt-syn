# System

Each type of {doc}`dynamical system <../../usage/problem_formulation/systems/index_systems>` has a dedicated  collection of routines to pose the Analysis and Synthesis problems. 


The routines are 
1. The `opt_system` algorithmic interconnection,
2. The `regulator` to build an internal model,
3. `lmi_analysis` and `lmi_synthesis` objects to pose the required LMIs.


The supported types of dynamical systems are:
```{toctree}
:maxdepth: 1
Linear Time Invariant <doc_lti>
Switched <doc_switched>
Periodic <doc_periodic>
Periodic-Orbit <doc_periodic_orbit>
``` 



All component routines inherit from the a `generic` interface.

## System
<!-- ```{eval-rst}
.. mat:autoclass :: system.generic
    :members:
``` -->

