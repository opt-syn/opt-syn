# Systems

Each type of {doc}`dynamical system <../../usage/problem_formulation/system/index_system>` has a dedicated  collection of routines to pose the analysis and synthesis problems. 


The routines are 
1. The `opt_system` algorithmic interconnection,
2. The `regulator` to build an internal model,
3. `lmi_analysis` and `lmi_synthesis` objects to pose the required LMIs.



The supported types of dynamical systems are
```{toctree}
:maxdepth: 1
Linear time-invariant <doc_lti>
Switched <doc_switched>
Periodic <doc_periodic>
Periodic-orbit <doc_periodic_orbit>
``` 


All component routines inherit from the `generic` interface.
## System (algorithmic interconnection)

```{eval-rst}
.. mat:autoclass :: system.generic.opt_system_interface   
    :members:
```

## Regulator

Both analysis and synthesis require a verification of the regulator equation.
```{eval-rst}
.. mat:autoclass :: system.generic.regulator_interface
    :members:
```


For analysis, the output of the `check_regulator()` function is contained in a separate class. If `check_regulator()` fails, then the algorithmic interconnection is not guaranteed to converge.

```{eval-rst}
.. mat:autoclass :: system.generic.reg_cl_out
    :members:
```

## LMI handler

These routines are shared by both analysis and synthesis.
```{eval-rst}
.. mat:autoclass :: system.generic.lmi_dispatch_interface
    :members:
```

## LMI analysis 

```{eval-rst}
.. mat:autoclass :: system.generic.lmi_analysis_interface
    :members:
    :show-inheritance:
```


## LMI synthesis 

```{eval-rst}
.. mat:autoclass :: system.generic.lmi_synthesis_interface
    :members:
    :show-inheritance:
```



## Dissipation container

The dissipation constraints are stored in `diss_data`. This container is used to construct the LMIs in analysis and synthesis.
```{eval-rst}
.. mat:autoclass :: manager.diss_data   
    :members:
```