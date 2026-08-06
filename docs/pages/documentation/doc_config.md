# Configuration

`opt_config` is the configuration file for the {doc}`Manager <doc_manager>` classes.

```{eval-rst}
.. mat:autoclass :: config.opt_config     
    :members:   
```

The specific subconfiguration options are:

## Numerical Tolerances

```{eval-rst}
.. mat:autoclass :: config.opt_config_tol 
    :members:   
```

## All programs (generic)

```{eval-rst}
.. mat:autoclass :: config.opt_config_gen    
    :members:   
```

## Analysis
```{eval-rst}
.. mat:autoclass :: config.opt_config_ana 
    :members:   
```

## Synthesis

```{eval-rst}
.. mat:autoclass :: config.opt_config_syn
    :members:   
```


## Bisection
The bisection options are used only when the Analysis or Synthesis problems are solved in Bisect or Alternating mode. These options are not required if the problem is solved only once, such as at a fixed linear rate $\rho$. 

```{eval-rst}
.. mat:autoclass :: config.bisect_opts        
   :members:
```