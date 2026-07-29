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
 <!-- so a specification is
always paired with a performance channel.
The properties themselves (linear convergence, input-to-state stability, $\ell_2$-gain)
are defined precisely on the {doc}`Performance Specifications <../../documentation/doc_specs>`
page. -->
Most specifications have an additional field `rho` discount rate $\rho >0$ as an argument. 
Choosing $\rho < 1$ imposes that the property holds at an exponential rate. 


Refer to {doc}`Performance Specifications <../../documentation/doc_specs>` for information about the specifications and their interfaces. 



<!-- ## Usage

A performance specification for Analysis or Synthesis is built from a cell of individual specifications.

```matlab
% linear convergence at rate rho (no performance channel needed)
rho = 0.92;
spec1 = spec_stability(rho);

% l2 gain on an existing performance channel with indices (iwp, izp):
GAIN = 10;
spec2 = spec_e2e(GAIN, iwp, izp);

%Impose only linear convergence
specs = {spec1};

%Impose both
specs = {spec1, spec2};

% l2 gain at 
spec2_rho = spec_e2e(GAIN, iwp, izp);
spec2.rho = rho;
specs = 
``` -->

<!-- A specification is built, optionally flagged as the optimization `target`, and passed to -->
<!-- an Analysis or Synthesis {doc}`manager <../../documentation/doc_manager>` together with the -->
<!-- IQC `order`:

```matlab
rho = 0.92;

% linear convergence at rate rho (no performance channel needed)
perf1 = spec_stability(rho);

% ...or an l2 gain on an existing performance channel with indices (iwp, izp):
perf2 = spec_e2e(GAIN, iwp, izp);
perf2.target = true;             % minimize the gain during the solve

% solve once,
sol = man_ana.solve_single(order, perf);

%or bisect for the best rate/gain
[sol_best, v_range] = man_ana.bisect(order, perf, bisect_opts());
``` -->

## Linear Convergence

Requests exponential stability of the iterates at rate $\rho$. For every
initial condition $x_0$ there is a fixed point $x^*(x_0)$ with

```{math}
\begin{align*}
\mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2
  \leq \gamma_0\, \rho^{k} \norm{x_0 - x^*(x_0)}_2, & & \forall k \in \N.
\end{align*}
```

Linear convergence holds if $\rho < 1$. 

This is the most common specification and needs no performance channel of its own (by default `iwp=[], izp=[]`).

It is specified with
```matlab
perf = spec_stability(rho);
```

:::{note}
If no specifications are supplied (`specs = []`), then the default specification used is `spec_stability(1)` with `target=true`.
:::


## Quadratic Performance

Several  specifications are special cases of general **quadratic supply-rate** conditions on
$(w_p, z_p)$. A quadratic supply rate condition with respect to matrices $(Q, S, R)$ is the existence of an $\epsilon > 0$ such that 

```{math}
\begin{align*}
\sum_{k=0}^{T}
 \rho^{-2k} \mat{c}{w_{p,k} \\ z_{p,k}}^\top
  \mat{cc}{Q & S \\ S^\top & R}
  \mat{c}{w_{p,k} \\ z_{p,k}} \; \preceq \; \sum_{k=0}^T -\epsilon \rho^{-2k} \norm{w_{p, k}}_2^2
\end{align*}
```
for all time horizons $T > 0$. 

Quadratic performance specifications are specified by 
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

If $\rho \in (0, 1)$, then &#8467;2 stability implies an Input-to-State Stability property 
<!-- This  which constitutes an input-to-state stability condition. There is no performance -->
<!-- output; only a bounded penalty on $w_p$ is imposed, giving -->
```{math}
\begin{align*}
\norm{x_k - x^*(x_0)}_2^2 \leq \gamma_x\, \rho^{2k} \norm{x_0 - x^*(x_0)}_2^2
  + \frac{\gamma \rho^2}{1-\rho^2} \max_{t \in 0, \ldots, k} \norm{w_{p, t}}_2^2, & & \forall k \in \N.
\end{align*}
```


Linear convergence is recovered from Input-to-State Stability when the
disturbance vanishes ($w_p= 0$). 


Reach for this when you need the iterates to stay
bounded and convergent under noise but do not need a specific gain; for a gain bound, use
$\ell_2$ Gain instead.

```matlab
perf = spec_l2(iwp);        % no performance output; optional bound: spec_l2(iwp, MU)
```

### &#8467;2 Gain

*Energy-to-energy gain.* Bounds the induced $\ell_2$-gain $\gamma$ from the performance
input to the performance output,

```{math}
\begin{align*}
\limsup_{T \to \infty}
  \frac{\sum_{k=0}^{T} \rho^{2k} \norm{z_{p,k}}_2^2}{   \sum_{k=0}^{T} \rho^{-2k} \norm{w_{p,k}}^2_2}
  \; < \; \gamma^2.
\end{align*}
```

For a linear system this, this is the $H_\infty$ gain. 


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

Setting both indices to zero requests plain passivity; positive indices request the
correspondingly stronger input- or output-strict passivity properties. The performance input and
output channels must have the same length.

```matlab
% ind_w = nu_w, ind_z = nu_z
perf = spec_passivity(ind_w, ind_z, iwp, izp);   
```

## Ergodic Convergence

Ergodic convergence arises in the optimization setting where all operators are subdifferentials ($F_i = \partial f_i$), and their operators classes are {class}`op_sml`, {class}`op_pcc`, or {class}`op_quad`.

An algorithm with no repeated operator evaluations satisfies ergodic convergence if there exists a $\gamma>0$ such that 
```{math}
\begin{align*}
\sum_{i=1}^s \left[f_i(z_i) - f_i(z_i^*(x_0)) - (w^*_i)^\top (z^*_i - z^*_i)\right] \leq \frac{\gamma}{k+1} \norm{x_0 - x^*(x_0)}_2^2 & & \forall k \in \N.
\end{align*}
```

The bracketed quantity is equal to 0 at optimality. The IQC-based formulation of ergodic convergence is based on the work of {footcite}`upadhyaya2025automated` (Section 4.1.2).


Ergodic convergence requires the introduction of new performance channels for the $(w^*_i)^\top z^*_i$ term. It is implemented as 
```matlab
[perf_erg, sys] = spec_ergodic(sys);
specs = {perf_erg};
```

Ergodic convergence is weaker than linear convergence. It can certify properties of convex optimization algorithms, whereas establishment of global linear convergence requires strong convexity. Ergodic convergence should only be used if all performance specifications have $\rho=1$. 






:::{warning}
In the current implementation, Ergodic convergence requires nonstrict feasibility of linear matrix inequalities. In numerical experiments, the minimal eigenvalue of a positive-definite-constrained block is $\approx (-10^{-12})$, which is not greater than or equal to  $0$. Future developments will hopefully patch this feasibility issue, in the meantime use with caution.
:::


## More Specifications

We plan to implement further performance criteria for both Analysis and Synthesis.