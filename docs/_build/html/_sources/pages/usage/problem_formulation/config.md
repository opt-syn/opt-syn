# Configuration

The configuration options are called by 
``` matlab
config = opt_config();
```

Configuration options include numerical tolerances and  recovery prereferences. These are detailed in {doc}`Configuration <../../documentation/doc_config>`.


## Restricting  Information Structures

The sparsity pattern of the controller $D_K$ in Synthesis can be imposed in the configuration options. This sparsity pattern can be used to enforce an {doc}`Information Structure <../../how_it_works/alg_properties>` on the generated algorithm.

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

By default, `config.syn.D_mask` will be an lower-triangular matrix with all ones.


An algorithm with nonzero upper-block-triangular  entries of `D_mask` can be Analyzed or Synthesized. However {{osyn}} cannot guarantee that the resulting algorithm will be well-posed, nor will it be able to {doc}`Simulate <../simulation>` trajectories of an algorithm execution.