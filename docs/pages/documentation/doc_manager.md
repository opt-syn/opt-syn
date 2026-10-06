# Manager


The {class}`manager` class is the highest level class in the {{osyn}} project. The analysis and synthesis problems are posed using the {class}`opt_analysis` and {class}`opt_synthesis` classes, respectively. The {class}`manager` classes are invoked in the {doc}`Problem formulation <../usage/problem_formulation/index_problem_formulation>` page, and their usage is explained in the {doc}`Solve <../usage/solve>` page.



Both the analysis and synthesis routines inherit from {class}`opt_manager_interface`, containing the common methods. 


The core user-facing methods are {meth}`solve_single`, {meth}`bisect`, and {meth}`alternate` for synthesis. 

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
## Common routines

```{eval-rst}
.. mat:autoclass :: manager.opt_manager_interface   
    :members:
```