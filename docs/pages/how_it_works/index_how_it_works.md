# How it Works

 
{{osyn}}'s analysis and synthesis routines are based on links between optimization and robust control. These underlying principles are explained in the following sections:

```{toctree}
:maxdepth: 1
Optimization Algorithms <opt_algs>
Properties of Algorithms <alg_properties>
Convergence Conditions <convergence_conditions>
IQC Analysis <op_uncertainty_iqc>
Synthesis and Networks <network_synthesis>
```

This overview is limited to static optimization problems with time-independent memory/stepsize rules. The {doc}`Problem Formulation <../usage/problem_formulation/index_problem_formulation>` section in {doc}`Usage<../usage/index_usage>` documents generalizations to this base construction, including 
- Time-varying optimization problems
- Time-varying dynamical systems
- Repeated operator calls
