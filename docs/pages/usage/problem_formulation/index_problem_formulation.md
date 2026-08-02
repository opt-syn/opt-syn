# Problem Setup 

Analysis and Synthesis procedures are both specified by the same three properties
1.  {doc}`System <system/index_system>` (from the last pages)
2. Configuration options
3. Performance specifications

The System stores details about the operator classes, networks, and algorithms used to solve the inclusion problem. 


Configuration defines options such as numerical tolerances.

The Performance Specifications define the convergence rate and considered  robustness criteria  in Analysis and Synthesis.


The System and Configuration are used to define the 
Analysis and Synthesis   {doc}`managers <../../documentation/doc_manager>`:
```matlab
sys = opt_system([arguments]);
config = opt_config();
man_ana = opt_analysis(sys, config);  
man_syn = opt_synthesis(sys, config); 
```

The managers and specifications are subsequently used to {doc}`Solve <../solve>` the Analysis and Synthesis problems with respect to the Performance Specifications.

```{toctree}
:maxdepth: 1
Performance Specifications <specs>
Configuration <config>
```