# Performance Specifications

A performance specification is the property you ask {{osyn}} to certify about an
algorithm. In Analysis it is the property that is verified for a given algorithm; in
Synthesis it is the property the returned algorithm is required to satisfy. Every
specification is a condition on the performance channel $(w_p, z_p)$ of the
{doc}`generalized plant <../../documentation/plants/doc_genplant>`, so a specification is
always paired with a performance channel.

The properties themselves (linear convergence, input-to-state stability, $\ell_2$-gain)
are defined precisely on the {doc}`Performance Specifications <../../documentation/doc_specs>`
page.

Most specifications accept a discount rate $\rho \in (0, 1]$. Choosing $\rho < 1$ requests
the property at an exponential rate: the specification is then evaluated on the
$\rho$-weighted signals, and every other specification on the same problem inherits the
same weighting.

## Usage

A specification is built, optionally flagged as the optimization `target`, and passed to
an Analysis or Synthesis {doc}`manager <../../documentation/doc_manager>` together with the
IQC `order`:

```matlab
rho = 0.92;

% linear convergence at rate rho (no performance channel needed)
perf = spec_stability(rho);

% ...or an l2 gain on an existing performance channel (iwp, izp):
perf = spec_e2e(GAIN, iwp, izp);
perf.target = true;             % minimize the gain during the solve

% solve once, or bisect for the best rate/gain
sol = man_ana.solve_single(order, perf);
% [sol_best, v_range] = man_ana.bisect(order, perf, bisect_opts());
```

Set `target = true` on the quantity you want optimized (the gain, or the rate under
bisection); leave it `false` to check a fixed bound. The performance-channel indices
`iwp`, `izp` come from the generalized plant definition.

## Linear Convergence

Requests exponential stability (linear convergence) of the iterates at rate $\rho$. For every
initial condition there is a fixed point $x^*(x_0)$ with

```{math}
\begin{align*}
\mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2
  \leq \gamma_0\, \rho^{k} \norm{x_0 - x^*(x_0)}_2, & & \forall k \in \N.
\end{align*}
```

This is the most common specification and needs no performance channel of its own. 
It is specified with

```matlab
perf = spec_stability(rho);
```

## Quadratic Performance

Some of the specifications are special cases of **quadratic supply-rate** conditions on
$(w_p, z_p)$, which imposes

```{math}
\begin{align*}
\sum_{k=0}^{T}
  \mat{c}{z_{p,k} \\ w_{p,k}}^\top
  \mat{cc}{Q & S \\ S^\top & R}
  \mat{c}{z_{p,k} \\ w_{p,k}} \; \preceq \; 0
\end{align*}
```

for a particular choice of the multiplier $\mat{cc}{Q & S \\ S^\top & R}$. 
It is specified with

```matlab
M = [Q, S; S', R];
perf = spec_quad(M, iwp, izp);
```

### &#8467;2 Gain

*Energy-to-energy gain.* Bounds the induced $\ell_2$-gain $\gamma$ from the performance
input to the performance output,

```{math}
\begin{align*}
\limsup_{T \to \infty}
  \frac{\sum_{k=0}^{T} \norm{z_{p,k}}_2}{\sum_{k=0}^{T} \norm{w_{p,k}}_2}
  \; < \; \gamma .
\end{align*}
```

For a linear system this is the $H_\infty$ gain. Use it to quantify how strongly a
disturbance at $w_p$ (for example, noise in the gradient evaluations) is amplified at a
tracked output $z_p$. Flag it as the `target` to minimize the gain.

```matlab
perf = spec_e2e(GAIN, iwp, izp);   % GAIN is the initial bound
perf.target = true;
```

### &#8467;2 Stability

*Stability and ISS.* Certifies that the performance input has a bounded effect on the
state, which constitutes an input-to-state stability condition. There is no performance
output; only a bounded penalty on $w_p$ is imposed, giving

```{math}
\begin{align*}
\norm{x_k - x^*(x_0)}_2^2 \leq \gamma_x\, \rho^k \norm{x_0 - x^*(x_0)}
  + \gamma_w \max_{t \in 0, \ldots, k} \norm{\delta w_t}_2^2, & & \forall k \in \N.
\end{align*}
```

With $\rho \in (0, 1]$ this establishes ISS, and it recovers linear convergence when the
disturbance vanishes ($\delta w_k = 0$). Reach for this when you need the iterates to stay
bounded and convergent under noise but do not need a specific gain; for a gain bound, use
$\ell_2$ Gain instead.

```matlab
perf = spec_l2(iwp);        % no performance output; optional bound: spec_l2(iwp, MU)
```

### Passivity

Imposes a passivity relation between the performance input and output, optionally with an
input passivity index $\nu_w$ and output passivity index $\nu_z$: for all horizons with
$x = 0$,

```{math}
\begin{align*}
\sum_{k=0}^{T} z_{p,k}^\top w_{p,k}
  \; \geq \;
  \sum_{k=0}^{T} \big( \nu_w \norm{w_{p,k}}_2^2 + \nu_z \norm{z_{p,k}}_2^2 \big).
\end{align*}
```

Setting both indices to zero requests plain passivity; positive indices request the
correspondingly stronger input- or output-strict passivity. The performance input and
output channels must have the same length.

```matlab
perf = spec_passivity(ind_w, ind_z, iwp, izp);   % ind_w = nu_w, ind_z = nu_z
```

## More Specifications

We plan to implement further performance criteria for both Analysis and Synthesis,
including peak-to-peak (generalized $H_\infty$) and $H_2$-type covariance amplification.