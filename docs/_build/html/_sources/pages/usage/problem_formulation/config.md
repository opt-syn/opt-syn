# Configuration

The configuration options are called by 
``` matlab
config = opt_config();
```

Configuration options include numerical tolerances and  recovery prereferences. These are detailed in {doc}`Configuration <../../documentation/doc_config>`.


## Restricting  Information Structures

The sparsity pattern of the controller $D_K$ in Synthesis can be imposed in the configuration options. This sparsity pattern can be used to enforce an [Information Structure](project:#information-structure) on the generated algorithm.

As an example, a two-operator splitting algorithm with a forward evaluation on $F_1$ and a backward evaluation on $F_2$ can be synthesized using the configuration options
```matlab
config.syn.D_mask = [0, 0; 1, 1]; %F2 can use information from F1
config.syn.D_mask = [0, 0; 0, 1]; %F2 cannot use information from F1.
```
 
By default, `config.syn.D_mask` will be an block-lower-triangular matrix with all ones (all operators admit known backward evaluations).