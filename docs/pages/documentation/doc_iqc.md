# Integral Quadratic Constraints

Input-output sequences of the operators satisfy a family of relations defined by {doc}`Integral Quadratic Constraints (IQC) <../how_it_works/op_uncertainty_iqc>`.  Analysis searches over the IQCs, while Synthesis returns a controller with fixed IQCs.



## IQCS for Analysis

The IQC structure for analysis of operators is a tuple $(\Psi, M, X, \ell)$, in which $\Psi$ is structured as 

```{math}
\begin{align}
        \Psi = \mat{cc}{\Psi_1 & 0 \\ 0 & \Psi_2}: & & \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p^\psi_{k} \\ q^\psi_{k}} =  \mat{cc|cc}{A_1 & 0 & B_1 & 0 \\
        0 & A_2 & 0 & B_2 \hl
        C_1 & 0 & D_1 & 0 \\
        0 & C_2 & 0 & D_2 
}  \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p_k \\ q_k}.
\end{align}
```

The analysis IQCs may be tall: the dimensions of $(p^\psi, q^\psi)$ may be larger than the dimensions of $(p, q)$ respectively.

## IQCS for Synthesis

IQC synthesis requires filters in which the dimensions of $(p^\psi, q^\psi)$ and $(p, q)$ are equal. The class of filters used are

```{math}
\begin{align}
        \Psi = \mat{cc}{\Psi_1 & \Psi_3 \\ 0 & \Psi_2}: & & \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p^\psi_{k} \\ q^\psi_{k}} =  \mat{cc|cc}{A_1 & 0 & B_1 & 0 \\
        0 & A_2 & 0 & B_2 \hl
        C_1 & C_3 & D_1 & D_3 \\
        0 & C_2 & 0 & D_2 
}  \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p_k \\ q_k}.
\end{align}
```

In IQC synthesis, it is required that $\Psi_1, \Psi_2, \Psi_3$ are all stable, and that $\Psi_2$ has a stable inverse. A tall IQC from analysis must be factorized into this square form first.

## IQC Routines

The iqcs routines are `iqc_loop_split` (analysis) and `iqc_loop_factored` (synthesis).

### Analysis Structure

```{eval-rst}
.. mat:autoclass :: iqc.iqc_loop_split   
    :members:
```

### Synthesis Structure



```{eval-rst}
.. mat:autoclass :: iqc.iqc_loop_factored 
    :members:
```

### Data Container

```{eval-rst}
.. mat:autoclass :: iqc.iqc_data_container 
    :members:
```