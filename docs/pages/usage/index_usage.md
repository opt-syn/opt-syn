# Usage


This page summarizes the usage of {{osyn}}.  Further detail on the {{osyn}} workflow is available in subsequent pages:
```{toctree}
:maxdepth: 1
Build the System <problem_formulation/system/index_system>
Simulate <simulation>
Problem Setup <problem_formulation/index_problem_formulation>
Solve and Validate <solve>
```


All code is written in object-oriented MATLAB. The {doc}`Documentation <../documentation/index_documentation>` page contains details about each object and function.

<!-- 
## Tasks

The three main tasks of {{osyn}} are Simulation, Analysis, and Synthesis.

Simulation solves an inclusion algorithm $0 \in \sum_{i=1}^{s} F_i(\beta^*)$ by iteratively executing an algorithm. Analysis certifies worst-case performance specifications of the algorithm for any $F_i$ in given classes of operators. Synthesis creates a controller such that the overall algorithm obeys the desired performance specifications.





All three tasks model the algorithm as a System ({class}`opt_system`). The System has three main parts: the operators, the network, and the controller.
```matlab
sys = opt_system(Operators, Network, Controller);
```

The operators are specified as an $s$-length cell array. These are specific operators in Simulation, or are operator classes in Analysis/Synthesis. The network and controller are both state-space dynamical systems. 

The network can be left empty `[]` if there are no network dynamics (standard optimization setting).

## Simulation
Simulation is conducted by the {class}`alg_sim` object. The first $T$ time steps of executing a $d$-dimensional inclusion algorithm is accomplished by
```matlab
simulator = alg_sim(sys, d);
sim_result = simulator.sim(T);
```

The output of simulation can then be plotted by {class}`alg_plotter`
```matlab
plt = alg_plotter(sim_result);
fig1 = plt.plot(desired signals); %creates MATLAB figure
```

## Analysis and Synthesis


The Analysis program certifies performance specifications by searching over valid Integral Quadratic Constraint relations satisfied by the operators $\{F_i\}$. Synthesis searches over controllers to d


Analysis and Synthesis are interfaced through manager objects
```matlab
man_ana = opt_analysis(sys);
man_syn  = opt_synthesis(sys);
```

The Analysis and Synthesis problems are specified by the  warm-starts, configuration options, and performance specifications. They are solved by  using the {meth}`solve_single`, {meth}`bisect`, or {meth}`alternate` commands.

```matlab
sol = man_ana.bisect(order, specs); %cell array of specifications
sol = man_syn.bisect(); %empty: default to linear convergence with rate 1
sol = man_syn.alternate(Niteration, order); %alternate between Synthesis and Analysis
```

<!-- Synthesis is -->


<!-- Analysis is additionally  -->
The order is a cell array of orders for the IQC program


The Analysis/Synthesis problem is solved by




Alternation between Analysis and Synthesis for 

 -->
