# Systems

Each type of {doc}`dynamical system <../../usage/problem_formulation/systems/index_systems>` has a dedicated  collection of routines to pose the Analysis and Synthesis problems. 


The routines are 
1. The `opt_system` algorithmic interconnection,
2. The `regulator` to build an internal model,
3. `lmi_analysis` and `lmi_synthesis` objects to pose the required LMIs.



The supported types of dynamical systems are
```{toctree}
:maxdepth: 1
Linear Time Invariant <doc_lti>
Switched <doc_switched>
Periodic <doc_periodic>
Periodic-Orbit <doc_periodic_orbit>
``` 


All component routines inherit from the `generic` interface.
## System (Algorithmic Interconnection)

```{eval-rst}
.. mat:autoclass :: system.generic.opt_system_interface   
    :members:
```

## Regulator

Both Analysis and Synthesis require a confirmation of the Regulator Equation.
```{eval-rst}
.. mat:autoclass :: system.generic.regulator_interface
    :members:
```


For Analysis, the output of the `check_regulator()` function is contained in a separate class. If `check_regulator()` fails, then the algorithmic interconnection is not guaranteed to converge.

```{eval-rst}
.. mat:autoclass :: system.generic.reg_cl_out
    :members:
```

## LMI Handler

These routines are shared by both Analysis and Synthesis.
```{eval-rst}
.. mat:autoclass :: system.generic.lmi_dispatch_interface
    :members:
```

## LMI Analysis 

```{eval-rst}
.. mat:autoclass :: system.generic.lmi_analysis_interface
    :members:
    :show-inheritance:
```


## LMI Synthesis 

```{eval-rst}
.. mat:autoclass :: system.generic.lmi_synthesis_interface
    :members:
    :show-inheritance:
```



## Dissipation Container

The dissipation constraints are stored in `diss_data`. This container is used to construct the LMIs in Analysis and Synthesis.
```{eval-rst}
.. mat:autoclass :: manager.diss_data   
    :members:
```