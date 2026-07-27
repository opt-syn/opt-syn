# Configuration

`opt_config` is the configuration file for the {doc}`Manager <doc_manager>` classes.

```{eval-rst}
.. mat:autoclass :: config.opt_config     
    :members:   
```

## Sub-configuration

The dedicated sub-configuration options are

```{eval-rst}
.. mat:autoclass :: config.opt_config_tol 
    :members:   
```

```{eval-rst}
.. mat:autoclass :: config.opt_config_gen    
    :members:   
```

```{eval-rst}
.. mat:autoclass :: config.opt_config_ana 
    :members:   
```

```{eval-rst}
.. mat:autoclass :: config.opt_config_syn
    :members:   
```


## Bisection
`bisect_opts` is used only when the Analysis or Synthesis problems are solved in Bisect or Alternating mode. `bisect_opts` is not required if the problem is solved only once, such as at a fixed linear rate $\rho$. 
```{eval-rst}
.. mat:autoclass :: config.bisect_opts        
   :members:
```