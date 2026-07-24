# Problem Setup 

Analysis and Synthesis procedures are both specified by the same three properties:
```{toctree}
:maxdepth: 1
System <system/index_system>
Configuration <config>
Performance Specifications <specs>
```

The System stores details about the operator classes, networks, and algorithms used to solve the inclusion problem. 


The Configuration options define options such as numerical tolerances.


The System and Configuration are used to define 
Analysis and Synthesis the  {doc}`managers <../../documentation/doc_manager>`:
```matlab
sys = opt_system([arguments]);
config = opt_config();
man_ana = opt_analysis(sys, config);
man_syn = opt_synthesis(sys, config);
```



The Performance Specifications define what convergence rate and robustness criteria are considered in algorithm Analysis and Synthesis.

 The managers and specifications are subsequently used to {doc}`Solve <../solve>` the Analysis and Synthesis problems with respect to the Performance Specifications.