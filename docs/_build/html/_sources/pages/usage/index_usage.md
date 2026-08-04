# Usage

These pages collect together the functionality of  {{osyn}}.  

The three main tasks are algorithm Simulation, Analysis, and Synthesis. 


Simulation solves an inclusion algorithm $0 \in \sum_{i=1}^{N_s} F_i(\beta^*)$ by iteratively executing a dynamical procedure (algorithm). Analysis certifies worst-case performance specifications of the algorithm for any $F_i$ in given classes of operators. Synthesis creates a controller that obeys the performance specifications.

<!-- All three tasks rely on the core construction of the System, or algorithmic interconnection.  -->

```{toctree}
:maxdepth: 1
:hidden:
Build the System <problem_formulation/system/index_system>
Simulate <simulation>
Problem Setup <problem_formulation/index_problem_formulation>
Solve and Validate <solve>
```


All code is written in object-oriented MATLAB. The {doc}`Documentation <../documentation/index_documentation>` page contains details about each object and function.


