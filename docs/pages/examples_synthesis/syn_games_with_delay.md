# Games with Delay

This example continues the Douglas-Rachford Game {doc}`simulation <../examples_simulation/sim_dr_game>` demonstration. 

In this two-operator example, the pseudogradient $F_1$ is $1.4785$-monotone and $0.1605$-cocoercive, and the normal cone $F_2 = \partial \mathbb{I}_{\norm{\cdot}_{\infty} \leq 10}$ is the subdifferential of a p.c.c. function. The operators are therefore described by {class}`op_gen` and {class}`op_pcc` respectively.


Three algorithms to solve this inclusion problem are synthesized. 

1. Backward evaluation of both $F_1$ and $F_2$,
2. Forward evaluation of $F_1$, backward evaluation of $F_2$,
3. A delay of one time step before and after evaluation of $F_1$.

Synthesis is performed for each case without a warm start (`iqc=[]`). The certified convergence rates are $\rho < 0.7163$ for Algorithm 1, $\rho < 0.8734$ for Algorithm 2, and $\rho < 0.9430$ for Algorithm 3.

All three are certifiably  convergent, but Algorithm 1 has the least worst-case convergence rate.
Each Synthesis is accompanied by an algorithm Simulation starting from an initial condition $x_0=0$, empirically demonstrating this speed of convergence

Figure  [1](#game-bw) plots a trajectory of the backward evaluation algorithm.
:::{figure} _static/game_sim_bw_dark.png
:align: center
:class: only-dark
:name: game-bw
*Figure 1:* Algorithm with both backward evaluations
:::

:::{figure} _static/game_sim_bw_dark.png
:align: center
:class: only-light
:name: game-bw
*Figure 1:* Algorithm with both backward evaluations
:::

Figure  [2](#game-fw) plots a trajectory of the forward and backward evaluation algorithm.
:::{figure} _static/game_sim_fw_dark.png
:align: center
:class: only-dark
:name: game-fw
*Figure 2:* Algorithm with one forward and one  backward evaluation
:::

:::{figure} _static/game_sim_fw_dark.png
:align: center
:class: only-light
:name: game-fw
*Figure 2:* Algorithm with one forward and one  backward evaluation
:::

Figure  [3](#game-delay) plots a trajectory of the algorithm with time delay.
:::{figure} _static/game_sim_delay_dark.png
:align: center
:class: only-dark
:name: game-delay
*Figure 3:* Algorithm with delay on $F_1$
:::

:::{figure} _static/game_sim_delay_dark.png
:align: center
:class: only-light
:name: game-delay
*Figure 3:* Algorithm with delay on $F_1$
:::


```{literalinclude} ../../../examples/synthesis/game_synth_sim.m
:caption: Code for Nash Equilibrium Seeking
:language: matlab
:linenos:  true
```

