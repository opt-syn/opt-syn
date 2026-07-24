# Manager


The `manager` class is the highest level class in the {{osyn}} project. The Analysis and Synthesis problems are posed using the `opt_analysis` and `opt_synthesis` classes, respectively. The `manager`classes are invoked in the {doc}`Problem Formulation <../usage/problem_formulation/index_problem_formulation>` page, and their usage is explained in the {doc}`Solve <../usage/solve>` page.



Both the analysis and synthesis routines inherit from `opt_manager_interface`, containing the common methods. 


The core user-facing methods are `solve_single`, `bisect`, and (for Synthesis) `alternate`. 

## Analysis
```{eval-rst}
.. mat:autoclass :: manager.opt_analysis   
    :members:
```

## Synthesis
```{eval-rst}
.. mat:autoclass :: manager.opt_synthesis   
    :members:    
```
## Common Routines

```{eval-rst}
.. mat:autoclass :: manager.opt_manager_interface   
    :members:
```