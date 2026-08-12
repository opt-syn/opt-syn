# Performance Specifications


Performance specifications describe constraints on the behavior of the algorithm. Analysis attempts to certify that an algorithm obeys the specification, and Synthesis tries to form an algorithm that obeys the specifications. 




Every specification is a condition on the performance channel $(w_p, z_p)$ of the
{doc}`generalized plant <../../documentation/plants/doc_genplant>`.

All specifications have the fields
:::{list-table} 
:header-rows: 1
*   - Field
    - Description    
*   - `iwp`
    - indices of performance input $w_p$
*   - `izp`
    - indices of performance output $z_p$
*   - `target`
    - should this specification be minimized
:::
Most specifications have an additional field `rho` discount rate $\rho >0$ as an argument. 
Choosing $\rho < 1$ imposes that the property holds at an exponential rate. 


Refer to the  {doc}`Performance Specification Documentation <../../documentation/doc_specs>` for information about the specifications and their interfaces. 


## Linear Convergence

This imposes exponential stability of the iterates at rate $\rho$. For every
initial condition $x_0$ there is a fixed point $(x^*(x_0), w^*(x_0), z^*(x_0))$ with

```{math}
\begin{align*}
\mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2
  \leq \gamma_0\, \rho^{k} \norm{x_0 - x^*(x_0)}_2 = 0, & & \forall k \in \N.
\end{align*}
```

Linear convergence holds if $\rho < 1$. 

This is the most common specification and needs no performance channel (by default `iwp=[], izp=[]`).

It is specified with
```matlab
perf = spec_stability(rho);
```

:::{note}

The default specification used if none are supplied (`specs = []`) is `{spec_stability(1)}` with `target=true`.
:::


Infinite-horizon exponential stability is that for $x_0$, 
```{math}
\begin{align*}
\lim_{k \rightarrow \infty}\rho^{-k}\mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2
\end{align*}
```

## Quadratic Performance

Several  specifications are special cases of general **quadratic supply-rate** conditions on
$(w_p, z_p)$. These conditions must be respected for all performance input sequence $(w_{p, k})_{k \in \N}$ with a finite $\rho$-weighted &#8467;2 norm  $(\sum_{k=0}^T -\epsilon \rho^{-2k} \norm{w_{p, k}}_2^2 < \infty)$.


 A quadratic supply rate condition with respect to matrices $(Q, S, R)$ is the existence of an $\epsilon > 0$ such that 

```{math}
\begin{align*}
\sum_{k=0}^{T}
 \rho^{-2k} \mat{c}{w_{p,k} \\ z_{p,k}}^\top
  \mat{cc}{Q & S \\ S^\top & R}
  \mat{c}{w_{p,k} \\ z_{p,k}} \; \leq  \; \sum_{k=0}^T -\epsilon \rho^{-2k} \norm{w_{p, k}}_2^2
\end{align*}
```
for all time horizons $T > 0$. 

Infinite-horizon quadratic performance is that 
```{math}
\begin{align*}
  \limsup_{T \rightarrow \infty} \sum_{k=0}^{T} \rho^{-2k} \mat{c}{w_{p,k} \\ z_{p,k}}^\top
  \mat{cc}{Q & S \\ S^\top & R}
  \mat{c}{w_{p,k} \\ z_{p,k}} \; +  \; \epsilon \rho^{-2k} \norm{w_{p, k}}_2^2 \leq 0.
\end{align*}
```

Quadratic performance specifications are imposed by 
```matlab
M = [Q, S; S', R];
perf = spec_quad(M, iwp, izp);
```

The Synthesis procedure requires that $Q = Q^\top$ and $R \prec 0$. 

The below specifications are all specific instances of quadratic performance:

### &#8467;2 Stability

*Stability and ISS.* Certifies that the performance input has a bounded effect on the
state.

```{math}
\begin{align*}
\sum_{k=0}^T \rho^{-2k} \norm{x_k - x^*(x_0)}_2^2 \leq \gamma\norm{x_0 - x^*(x_0)}_2^2
  + \gamma \sum_{k=0}^T \rho^{-2k} \norm{w_{p, k}}_2^2, & & \forall k \in \N.
\end{align*}
```

If $\rho \in (0, 1)$, then &#8467;2 stability implies an Input-to-State Stability property {footcite}`sontag1989smooth`, {footcite}`schwenkel2026multi` 
```{math}
\begin{align*}
\norm{x_k - x^*(x_0)}_2^2 \leq \gamma\, \rho^{2k} \norm{x_0 - x^*(x_0)}_2^2
  + \frac{\gamma \rho^2}{1-\rho^2} \max_{t \in 0, \ldots, k} \norm{w_{p, t}}_2^2, & & \forall k \in \N.
\end{align*}
```


Linear convergence is recovered from Input-to-State Stability when the
disturbance vanishes ($w_p= 0$). 


Use this when you need the iterates to stay
bounded and convergent under noise but do not need a specific gain; for a gain bound, use
$\ell_2$ gain instead.

```matlab
perf = spec_l2(iwp);        % no performance output; optional bound: spec_l2(iwp, MU)
```

Infinite-horizon  &#8467;2 stability is 
```{math}
\begin{align*}
\limsup_{T \to \infty}
  \frac{\sum_{k=0}^{T} \rho^{-2k} \norm{x_k - x^*(0)}_2^2}{   \sum_{k=0}^{T} \rho^{-2k} \norm{w_{p,k}}^2_2}
  \; < \; \gamma^2.
\end{align*}
```

### &#8467;2 Gain

*Energy-to-energy gain.* Bounds the induced &#8467;2 gain $\gamma$ from the performance
input to the performance output,

```{math}
\begin{align*}
\limsup_{T \to \infty}
  \frac{\sum_{k=0}^{T} \rho^{2k} \norm{z_{p,k}}_2^2}{   \sum_{k=0}^{T} \rho^{-2k} \norm{w_{p,k}}^2_2}
  \; < \; \gamma^2.
\end{align*}
```

For a linear system this, this is the $H_\infty$ gain under a $\rho$-weighting. The &#8467;2 gain is an infinite-horizon penalty.


Use it to quantify how strongly a
disturbance at $w_p$ (for example, noise in the gradient evaluations) is amplified at a
tracked output $z_p$.

```matlab
perf = spec_e2e(GAIN, iwp, izp);   % GAIN is the initial bound
perf.target = false; %enforce the gain bound

perf.target = true; %minimize the gain 
```

### Passivity

Imposes a passivity relation between the performance input and output. Passivity may optionally include an 
input passivity index $\nu_w$ and an output passivity index $\nu_z$: 

Passivity is obeyed if for all time horizons $T$ with  with $x_0 = 0, x^*(x_0) = 0$, it holds that
```{math}
\begin{align*}
\sum_{k=0}^{T} z_{p,k}^\top w_{p,k}
  \; \geq \;
  \sum_{k=0}^{T} \big( \nu_w \norm{w_{p,k}}_2^2 + \nu_z \norm{z_{p,k}}_2^2 \big).
\end{align*}
```

Setting both indices to zero requests standard passivity; positive indices request the
correspondingly stronger input- or output-strict passivity properties. The performance input and
output channels must have the same length.

```matlab
% ind_w = nu_w, ind_z = nu_z
perf = spec_passivity(ind_w, ind_z, iwp, izp);   
```

## Stochastic Sensitivity

Stochastic sensitivity imposes a mean-square boundedness requirement on the performance output {footcite}`van2021speed`. The performance input sequence $\{w_{p, k}\}$ is a sequence of i.i.d. random variables. These inputs are zero-mean and bounded: there exists a known $\Omega \succ 0$ such that  $\E[w_{p, k}] = 0$ and $\E[w_{p, k} w_{p, k}^\top] \preceq \Omega$ at all $k \in \N$. The algorithm achieves stochastic sensitivity with gain $\gamma \geq 0$ if for all initial conditions $x_0$ and performance inputs $w_p$, it holds that 
```{math}
 \limsup_{T \rightarrow \infty} \frac{1}{T} \sum_{k=0}^{T}
        \mathbb{E}[|| z_{p,k}||^2_2] \leq \gamma^2 \text{Trace}(\Omega).
```

Stochastic sensitivity is invoked by the command
```matlab
perf = spec_h2(GAIN, Omega, iwp, izp);   %GAIN = gamma
```

Stochastic sensitivity is only certified if  $\rho = 1$ and the oracle input $z$ is independent of $w_p$. In contrast, the $\ell_2$ gain is usable if these conditions are violated.


 
## Ergodic Convergence

Ergodic convergence arises in the optimization setting where all operators are subdifferentials ($F_i = \partial f_i$), and their operators classes are {class}`op_sml`, {class}`op_pcc`, or {class}`op_quad`.

An algorithm with no repeated operator evaluations satisfies ergodic convergence if there exists a $\gamma>0$ such that 
```{math}
\begin{align*}
\sum_{i=1}^s \left[f_i(z^i_k) - f_i(z^{*,i}(x_0))\right] - (w^*)^\top (z_k - z^{*}(x_0)) \leq \frac{\gamma}{k+1} \norm{x_0 - x^*(x_0)}_2^2 & & \forall k \in \N.
\end{align*}
```
The left-hand side  is a duality gap, and is equal to 0 at optimality. 
This  formulation of ergodic convergence in duality gap originates from  {footcite}`upadhyaya2025automated` (Section 4.1.2).


Ergodic convergence is called by 
```matlab
spec = spec_ergodic();
```

Ergodic convergence is weaker than linear convergence. It can certify properties of convex optimization algorithms, whereas establishment of global linear convergence requires strong convexity. Ergodic convergence should only be used if all performance specifications have $\rho=1$. 

:::{warning}
In the current implementation, Ergodic convergence requires nonstrict feasibility of linear matrix inequalities. In numerical experiments, the maximal eigenvalue of a negative-semidefinite-constrained block is $\approx 10^{-12}$, which is not less than or equal to  $0$. Future developments will try to patch this feasibility issue, in the meantime use with caution.
:::




## Performance for Time-Varying Dynamical Systems

The specific performance constraints imposed by `specs` may vary for {doc}`systems with time-variations <system/dynamics>`.

All specifications on this page are used as presented for LTI, periodic, and periodic-orbit systems. 

For switched systems, a performance specification in `specs` imposes a  worst-case bound over all possible switching sequences.


## More Specifications

We plan to implement further performance criteria for both Analysis and Synthesis.