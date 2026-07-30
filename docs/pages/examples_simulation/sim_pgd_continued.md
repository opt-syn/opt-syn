# Projected Gradient Descent (Continued)


We continue simulation of the Projected Gradient Descent algorithm, see {doc}`<../get_started>` for the first part.

PGD is described by 

```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma I & -\gamma I \hl I &0 & 0 \\
 I & -\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial f(z_k^1) \\  \partial I_{\mathcal{Z}}(z_k^2)}
\end{align*},
```