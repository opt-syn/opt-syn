# Solve and Validate

:::{caution} 
Under Construction
:::
## Solve

### Single Solve

### Bisection

### Alternation

## Validation

The `sol` structure contains information about the solution of analysis/synthesis. The solution is feasible if the following conditions are met

| Name   |  Description  | Valid Condition |
|----| ---- | ----- | 
| `STATUS` | Feasibility of problem | `STATUS`=0 if feasible, `STATUS`$\neq$0  if infeasible |
| `dia` | Constraint violation | `dia`<0 if strictly feasible, `dia`=0 if marginally feasible, `dia` > 0 if infeasible |
| `gain` | Input passivity index and $H_\infty$ gain | Feasible if `gain(1)` < 0 and `gain(2)` < 1 |

If all of the above conditions are met, then linear convergence is established if and only if `sol.rho` < 1. A finite `sol.rho` > 1 establishes a bounded rate of divergence. No conclusions can be drawn about linear convergence if `sol.rho` = 1. 