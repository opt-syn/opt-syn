# Time-Varying Optimal Solutions


A trajectory $\{\beta^*_k\}$ is a critical path {footcite}`bianchin2026internal` of a time-varying inclusion problem if
```{math}
\begin{align}
0 \in \sum_{i=1}^{N_s} F_{ik}(\beta^*_k) & & \forall k \in \N.
\end{align}
```

Time-variation in the operators $F_k$ may arise from tracking a moving target.


The time-varying inclusion problem may be expressed as the existence of a pair $(\beta^*, w^*)$ with 
```{math}
\begin{align}
0 \in \sum_{i=1}^{N_s} w^i_k = 0, \qquad  w^i \in  F_{ik}(\beta^*_k) & & \forall k \in \N.
\end{align}
```

{{osyn}} supports time-variation in $\beta^*$ according to linear dynamics. Time-variation in  $w^*$ is not yet supported. This restriction is equivalent to the existence of time-independent operators $F_\bullet$ such that 
```{math}
\begin{align}
 F_k(\beta) = F_{\bullet i}(\beta - \beta^*_k) & & \forall k \in \N.
\end{align}
```


The linear system (signal generator) governing the path $\beta^*$ is described by 
```{math}
\begin{align}
 \mat{c}{\eta_{k+1} \\ \beta^*_k} = \mat{c}{S_\beta \\ R_\beta} \eta_{k}
\end{align}
```

Tracking of an the optimal solution is accomplished by setting the `tracking` field in `opt_system` to a struct with fields (`Sbeta`, `Rbeta`). 

Time-variation of `w^*` is not yet supported.

:::{seealso}
The {doc}`Tracking <../../../examples_simulation/sim_tracking>` example executes an algorithm with an oscillating optimal trajectory, with provided code.
:::
