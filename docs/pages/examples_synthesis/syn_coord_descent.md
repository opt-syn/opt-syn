# Coordinate Descent

This examples Synthesizes cyclic coordinate descent algorithms, continuing from the previous {doc}`simulation <../examples_simulation/sim_coord_descent>` and {doc}`Analysis <../examples_analysis/ana_coord_descent>` examples.




Two types of cyclic coordinate descent algorithms are the considered in this demonstration:
1. Full Knowledge: the entire vector  $\nabla f(\beta_k)$ is known,
2. Partial Knowledge: only the actively updated portion $[\nabla f(\beta_k)]^i$.

The reference coordinate descent algorithm in the prior examples is 
```{math}
\begin{align}
    \beta^{i(k)}_{k+1} = \beta^{i(k)}_{k} -\gamma \   [\nabla f(\beta_k)]^i.
\end{align}
```
This reference algorithm has Partial Knowledge.

Synthesis is performed to design a $c$-block cyclic coordinate descent algorithm with respect to a function $f \in S_{1, 2}$. Figure [1](#coord-6) plots the upper bounds on the convergence rate $\rho$ as $c$ increases. Algorithms with Partial Knowledge (orange, top) have a higher worst-case convergence rate when compared against algorithms with Full Knowledge (blue, bottom) at $c>1$. The two synthesized algorithms are equal in the trivial case of $c=1$.



:::{figure} _static/coord_synth_1_2_dark.png
:align: center
:class: only-dark
:name: coord-6
*Figure 1:* Convergence rate of coordinate descent algorithm v.s. number of blocks
:::

:::{figure} _static/coord_synth_1_2_light.png
:align: center
:class: only-light
:name: coord-6
*Figure 1:* Convergence rate of coordinate descent algorithm v.s. number of blocks.
:::


```{literalinclude} ../../../examples/synthesis/syn_coord_descent.m
:linenos: true
:caption: Code for coordinate descent synthesis at $m=1, L=2, c=6$
:language: matlab
:lines: 1-26
```


<!-- :::{seealso}
{doc}`Analysis <../examples_analysis/ana_coord_descent>` of the coordinate-descent scheme.
::: -->
