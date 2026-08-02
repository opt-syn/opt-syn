# Solve and Validate



The  Analysis and Synthesis problems can be solved once  the appropriate {doc}`managers <../documentation/doc_manager>` are created:

```matlab
sys = opt_system([arguments]);
config = opt_config();
man_ana = opt_analysis(sys, config);  
man_syn = opt_synthesis(sys, config); 
```


Three solution modes are available: 
1. [Single Solution](#single-solution)
2. [Bisection](#bisection)
3. [Alternation](#alternation)

(#single-solution)=
### Single Solution


The  {meth}`solve_single` routine performs Analysis or Synthesis at given  specifications. The output of {meth}`solve_single` is a {doc}`solution <../documentation/doc_solutions>` structure {class}`opt_solution`. The fields of the solution structure are explained further  below in [Validation](#validation).

### Analysis 
The Analysis single solution routine is 
```matlab
sol_ana_single = man_ana.solve_single(order, specs);
```

The `order` input is an $s$-length cell array. Each entry  of the `order` cell array is one or two  nonnegative integers:
```{list-table}
:header-rows: 1
* - Order Length
  - Operator Classes     
* - 1
  - {class}`op_gen`, {class}`op_sml_causal`,  {class}`op_quad_causal`
* - 2
  - {class}`op_sml`, {class}`op_pcc`,  {class}`op_quad`
```

The order defines length of the filters in the {doc}`IQCs <../documentation/doc_iqc>`. A higher order requries more computation, but can lead to improved bounds. An example order declaration is


```matlab
op1 = op_gen([arguments]);
op2 = op_sml([arguments]);
op3 = op_quad([arguments]);
order = {1, [2, 0], [1, 1]};
```

The two-input order routines default to zero if only one input is passed. As an example for {class}`op_sml`, order `2` is equivalent to `[2, 0]`.

The `specs` is a cell array of {doc}`specifications <problem_formulation/specs>`. 


:::{note}
Analysis is performed with respect to a common certificate among all performance specifications in `specs` {footcite}`scherer1997multiobjective`. Requiring this common certificate may yield conservatism as compared to doing Analysis separately for each specification in `specs`. 

The common certificate is required if [Alternation](#alternation) between Analysis and Synthesis is performed.
:::

Synthesis is invoked by 
```matlab
sol_syn_single = man_syn.solve_single(iqc_init, specs);
```

The `iqc_init` input is an $s$-length cell array of {doc}`IQCs <../documentation/doc_iqc>` to warm-start the Synthesis process. A default input of `iqc=[]` will compute valid warm-starts using each {doc}`operator class's <../documentation/operators/doc_operators>` 
{meth}`create_iqc_identity` routine.


The objective in Analysis or Synthesis is set using the `target` field in the  {doc}`operator class's <../documentation/doc_specs>`. 
If `specs{j}.target = true`, then the specification in $i$ is minimized. At most one specification $j$ can have `specs{j}.target = true`.
If `specs{j}.target=false` for all specifications $j$, then a feasibility problem is solved. 


:::{caution}
Different specifications can have different rate $\rho$ in Synthesis if and only if 
1. the order for each {class}`op_gen` is `0`
2. the order for each {class}`op_sml`, {class}`op_pcc`, {class}`op_quad` is `[nonnegative, 0]`. 

If these conditions fail, then all specifications
- must have the same rate $\rho$ (`opt_config.gen.same_rho = true`), 
- must hold only in an infinite-horizon sense (in the current implementation).

The {class}`opt_manager` classes will check these conditions, and will override `same_rho` to true if not already set. This issue is due to the presence of noncausal IQC filters.
 <!-- See {doc}`<../documentation/doc_iqc>` for more detail about this issue. -->
:::

(#bisection)=
### Bisection

The routine {meth}`bisect` performs bisection on a single specification. Parameters of bisection are set using the `bisect_opts` {doc}`configuration <../documentation/doc_config>`. The index of the specification in the `spec` to minimize using bisection is set via `bisect_opts.spec_ind`, with a default of index of 1.

The outputs of  {meth}`bisect` are the {class}`opt_solution` structure and the two-element array `v_range`. The entries of `v_range` are the lower and upper bounds of the bisected parameter.

As an example, the Analysis and Synthesis routines for minimizing  the linear convergence rate $\rho$ subject to a $\rho$-weighted  $\ell_2$ gain bound of 100 is 
```matlab
spec_rho = spec_stability();
spec_l2 = spec_e2e(100);
specs = {spec_rho, spec_l2};
b_opts = bisect_opts;
b_opts.spec_ind = 1;

[sol_ana_bisect, v_ana] = man_ana.bisect(order, specs, b_opts);
[sol_syn_bisect, v_syn] = man_syn.bisect(iqc_init, specs, b_opts);
```
(#alternation)=
### Alternation

The {meth}`alternate` routine switches between Synthesis and Analysis. It solves a Synthesis problem with fixed IQCs to find a controller, and then solves Analysis with the fixed controller to find IQCs. The Analysis and Synthesis problems may include inner bisection steps.

An alternation routine with 3 Synthesis/Analysis steps and inner bisection is performed by 
```matlab
Niter = 3; 
[sol_syn_alternate, v_history] = man_syn.alternate(Niter, iqc_init, order, specs, b_opts);
```

The `v_history` output is a cell array with 2 rows and  Niter columns. Each entroy of the cell array stores the  lower and upper parameter bounds from bisection. The top row are the Synthesis bounds, and the bottom row are the Analysis bounds.


##  Validate

The {class}`opt_solution` structure contains information about the solution of analysis/synthesis. The solution is feasible if the following conditions are met:
| Name   |  Description  | Valid Condition |
|----| ---- | ----- | 
| `STATUS` | Feasibility of problem | `STATUS`$=0$ if feasible, `STATUS`$\neq 0$  if infeasible |
| `dia` | Constraint violation | `dia`$<0$ if strictly feasible, `dia`$=0$ if marginally feasible, `dia` $> 0$ if infeasible |
| `gain` | Input passivity index and $H_\infty$ gain | Feasible if `gain(1)` $< 0$ and `gain(2)` $< 1$ |
| `rho` | Convergence rate |  Linearly convergent if $\rho < 1$ |
| `regcl` | Closed-loop regulator equation | Nonempty struct with fields (`S`, `R`, `Pi`, `Gam`, `Phi`, `Th`)| 


If the algorithm has a block-lower-triangular information structure, and  Analysis or Synthesis is successful at some $\rho > 0$, then the algorithm is also well-posed.


Other attributes of {class}`opt_solution` include
:::{list-table} 
:header-rows: 1
*   - Field
    - Description    
*   - `vars`
    - Design variables in the Analysis/Synthesis problem
*   - `sys`
    - System that solves the inclusion problem
*   - `lmi_out`
    - Solution information from the numerical solver (LMILAB)
*   - `cert`
    - Analysis-and-Synthesis-specific certificates
:::
 
 