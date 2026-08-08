# Sequential Games



The four-player game is continued from the  {doc}`simulation <../examples_simulation/sim_dr_game>` example. Its  pseudogradient  map $F_1$ is affine, $1.4785$-monotone and $0.1605$-cocoercive. 

This example performs Synthesis of a Nash Equilibrium seeking algorithm in two scenarios:
1. Simultaneous: each agent updates its strategy $z^i_k$ at every time step $k$
2. Sequential: the agents take turns updating their strategies, $z^i_k$ changes only once per four time steps.

Both algorithms permit only explicit evaluations of the pseudogradient map $F_1$. 

The Simultaneous  setting involves an LTI system, and returns an algorithm with convergence rate $\rho \leq 0.8733$. The Sequential setting is modeled as a periodic-orbit system, and yields an algorithm with convergence rate $\rho \leq 0.9667$.


Figure  [1](#game-simul) plots a trajectory of the simultaneous game.
:::{figure} _static/syn_nash_uncons_simul_dark.png
:align: center
:class: only-dark
:name: game-simul
*Figure 1:* Simultaneous play of the game
:::

:::{figure} _static/syn_nash_uncons_simul_light.png
:align: center
:class: only-light
:name: game-simul
*Figure 1:* Simultaneous play of the game
:::

Figure  [2](#game-coord) plots a trajectory of the sequential game.
:::{figure} _static/syn_nash_uncons_coord_dark.png
:align: center
:class: only-dark
:name: game-coord
*Figure 2:* Sequential play of the game
:::

:::{figure} _static/syn_nash_uncons_coord_light.png
:align: center
:class: only-light
:name: game-fw
*Figure 2:* Sequential play of the game
:::

In the designed sequential algorithm, each agent $i$ has full knowledge of the entire pseudogradient vector $w_k$. The partial knowledge setting  involves  allowing only the actively playing agent $i$ access to $w^i_k$. Synthesis in this partial knowledge setting returns an infeasible controller $(\rho \geq 2)$.

```{literalinclude} ../../../examples/synthesis/syn_nash_eq_coord.m
:caption: Code for sequential Nash Equilibrium Seeking
:language: matlab
:linenos:  true
:lines: 1-29
```

