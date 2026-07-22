# Problem Setup 

Analysis and Synthesis problems are both defined by the same three properties:
1. The system (operator classes, network dynamics)
2. Performance specifications,
3. Configuration options.



The system is specified by five attributes
``` matlab
sys = opt_system(Operators, Network, Controller, Bind, Tracking);
```

The configuration options are called by 
``` matlab
config = opt_config();
```

Configuration options include numerical tolerances, verbosity levels



```{toctree}
:maxdepth: 1
:hidden:
Operators <operators>
Dynamical Systems <systems/index_systems>
Bind <bind>
Tracking <tracking>
Performance Specifications <specs>
```

