<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/opt-syn/opt-syn/blob/main/docs/opt_syn_dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/opt-syn/opt-syn/blob/main/docs/opt_syn_light.svg">
  <img alt="opt-syn: Analysis and Synthesis of First-Order Algorithms">
</picture>

# opt-syn

Analysis and Synthesis of First-Order Algorithms (MATLAB)


<span style="color:#5D93BF">opt</span>-syn: <span style="color:#5D93BF">optimization</span> (algorithm) synthesis.

## Overview


<span style="color:#5D93BF">opt</span>-syn  analyzes and synthesizes first-order algorithms by using methods from robust control theory. 


Analysis tasks include establishment of worst-case linear convergence rates and  gains for error amplification. Synthesis tries to design an optimization algorithm satisfying these desired constraints.  These Analysis and Synthesis tasks may be performed for algorithms arising in dynamic environments,  including cases  with constant or time-varying delays, channel memory, and cross-talk.



Both the Analysis and Synthesis problems are posed as convex problems with  Linear Matrix Inequality constraints, and are solved using [LMILab](https://www.mathworks.com/help/robust/ug/introduction.html).  
