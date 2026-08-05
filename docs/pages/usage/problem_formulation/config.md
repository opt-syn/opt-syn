# Configuration

The configuration options are called by 
``` matlab
config = opt_config();
```

Configuration options include numerical tolerances and  recovery prereferences. These are detailed in the  {doc}`Configuration Documentation <../../documentation/doc_config>`.


## Restricting  Information Structures

The primary user-facing configuration option in algorithm synthesis is specification of the sparsity pattern of the controller matrix $D_K$. This sparsity pattern can be used to enforce an {doc}`Information Structure <../../how_it_works/alg_properties>` on the generated algorithm.

As an example, a two-operator splitting algorithm with $\Dcl_{11}=0$ (explicit evaluation of $F_1$) and $\Dcl_{22}\neq 0$ (implicit evaluation of $F_2$) can be synthesized using the configuration options
```matlab
%F2 can use information from F1
config.syn.D_mask = [0, 0; 
                     1, 1]; 

%F2 cannot use information from F1
config.syn.D_mask = [0, 0; 
                     0, 1]; 
```

Algorithms with both explicit evaluations ($\Dcl_{11}, \Dcl_{22} =  0$) can be synthesized using 
 ```matlab
 %F2 can use information from F1
config.syn.D_mask = [0, 0; 
                      1, 0];

%F2 cannot use information from F1 
config.syn.D_mask = [0, 0; 
                      0, 0]; 
```

By default, `config.syn.D_mask` will be an lower-triangular matrix with all ones, permitting all implicit evaluations.


An algorithm with nonzero upper-block-triangular  entries of `D_mask` can be Analyzed or Synthesized. However {{osyn}} cannot guarantee that the resulting algorithm will be well-posed, nor will it be able to {doc}`Simulate <../simulation>` trajectories of an algorithm execution.


## Simplified Synthesis Programs

Special structures of the IQC synthesis programs allow for simplification of the LMI programs. 
Supported simplification methods v.s. dynamical system types are
:::{list-table}
:header-rows: 1
:stub-columns: 1
* - LTI 
  - Periodic-Orbit
  - Periodic
  - Switched
* - Matrix Elimination
  - [x]
  - [x]
  - []
  - []
* - Reduced-Order
  - [x]
  - [x]
  - []
  - []
:::

All simplifications  are enabled by default. They may respectively be disabled by 
```matlab
config.syn.elimination = false;
config.syn.reduced_order = false;
```

### Matrix Elimination

If only one performance requirement is present, then the Matrix Elimination Lemma {footcite}`gahinet1994linear` may be used to remove the some or all of controller variables from the algorithm design problem. A single LMI constraint with the controller variables is replaced by two smaller LMI constraints lacking these variables. Both constraint sets have the same feasibility region.

The Synthesis programs with and without Elimination use the Transformation approach from {footcite}`scherer1997multiobjective` to design controllers.

### Reduced-Order Control

Reduced-Order Control uses the internal model structure of the controllers to lower the number of states in the generated algorithm. This reduced-order control is based on the formulation of {footcite}`korouglu2009generalized`. Reduced-order control also allows for the search over solutions of the Regulator Equations.