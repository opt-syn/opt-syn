# Nash Equilibrium Seeking

The Douglas-Rachford algorithm is a procedure for solving a two-operator inclusion problem {footcite}`douglas1956numerical`.
It is characterized by parameters $\gamma, \lambda > 0$, and can be described as the interconnection
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2)}.
\end{align*}
```


This example involves finding a variational Nash Equilibrium of a Linear-Quadratic game. 
Four agents are interacting in a noncooperative manner. 
The individual agent decisions are each described in a vector $\beta_v \in \R^{5}$ for $v \in \{1, \ldots, 4\}$. These decisions  are concatenated into a vector $\beta \in \R^{20}$. 

Each agent is trying to minimize their own quadratic cost 
```{math}
\begin{align*}
f_v(\beta) &= \frac{1}{2} \beta^\top Q_v \beta + b_v^\top \beta + e_v, & & \forall v \in \{1, \ldots, 4\}.
\end{align*}
```


The pseudogradient operator of this game is
```{math}
\begin{align*}
F_1 &= \mat{c}{\partial_{\beta_1} f_1 \\ \partial_{\beta_2} f_2 \\ \partial_{\beta_3} f_3 \\ \partial_{\beta_4} f_4}.
\end{align*}
```

The agent strategies are restricted such that $\norm{\beta}_{\infty} \leq 10$. 
The point $\beta^*$ is a variational Nash Equilibrium if it satisfies the variational inequality
```{math}
\langle F_1(\beta^*), \beta - \beta^* \rangle \geq 0, \qquad \forall \beta \ \text{with} \norm{\beta}_{\infty} \leq 10.
```
Under this condition, no agent $v$ will unilaterally change their decision $\beta_v$ to decrease their individual objective. The variational inequality can be converted to an inclusion problem with $F_2 = \partial \mathbb{I}_{\norm{\cdot}_\infty \leq 10}$

 $$0 \in F_1(\beta^*) + F_2(\beta^*).$$ 

In this example, the affine pseudogradient operator $F_1$ has strong monotonicity constant $\mu = 1.4785$ and cocoercivity constant $\beta = 0.1605$. This game therefore has a unique variational Nash Equilibrium $\beta^*$. 


We apply the Douglas Rachford algorithm with parameters $\gamma = 1, \lambda =1$ to find $\beta^*$. Figure [1](#dr-game) plots algorithm trajectories starting from the initial condition $x_0 = 0$. The bottom-left plot displays the payoffs $\{f_v(\beta_k)\}_{v =1}^4$ of the individual agents.

:::{figure} _static/dr_game_dark.png
:align: center
:class: only-dark
:name: dr-game
*Figure 1:* Convergence to the variational Nash Equlibrium
:::

:::{figure} _static/dr_game_light.png
:align: center
:class: only-light
:name: dr-game
*Figure 1:* Convergence to the variational Nash Equlibrium
:::


```{literalinclude} ../../../examples/simulation/dr_example/sim_dr_LQ_game.m
:linenos: true
:caption: Code for Douglas-Rachford based variational Nash Equilibrium seeking.
:language: matlab
```

