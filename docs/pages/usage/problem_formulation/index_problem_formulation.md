# Problem Setup 

Analysis and Synthesis problems are both defined by the same three properties:
1. System
2. Performance Specifications,
3. Configuration Options.


```{toctree}
:maxdepth: 1
:hidden:
Operators <operators>
Dynamical Systems <dynamical_systems/index_systems>
Bind <bind>
Tracking <tracking>
Performance Specifications <specs>
Configuration <config>
```


The System is specified by five attributes
``` matlab
sys = opt_system(Operators, Network, Controller, Bind, Tracking);
```

These attributes are explained in subsequent pages.


