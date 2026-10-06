# Problem setup 

Analysis and synthesis procedures are both specified by the same three properties
1.  {doc}`System <system/index_system>` (from the last pages)
2. Configuration options
3. Performance specifications

The system stores details about the operator classes, networks, and algorithms used to solve the inclusion problem. 


Configuration defines options such as numerical tolerances.

The performance specifications define the convergence rate and considered  robustness criteria  in analysis and synthesis.


The system and configuration are used to define the analysis and synthesis   {doc}`managers <../../documentation/doc_manager>`:
```matlab
sys = opt_system([arguments]);
config = opt_config();
man_ana = opt_analysis(sys, config);  
man_syn = opt_synthesis(sys, config); 
```

The managers and specifications are subsequently used to {doc}`Solve <../solve>` the analysis and synthesis problems with respect to the performance specifications.



```{toctree}
:maxdepth: 1
Performance specifications <specs>
Configuration <config>
```

:::{tip}
Omission of the `config` argument leads to use of the default configuration options in {class}`opt_config`:
```matlab
man_syn = opt_synthesis(sys); 
```
:::