# Generalized Plant


{#genplant}
The input-output signals of the plant are
:::{list-table} Input-Output Signals
:widths: 2 12 2 12
:header-rows: 1
*   - 
    - Plant Input
    -
    - Plant Output
*   - $z$
    - input to operators
    - $w$ 
    - output from operators
*   - $z_p$
    - performance output
    - $w_p$
    - Performance Input
*   - $y$
    - output to controller
    - $u$
    - input from controller
:::

<!-- *   - 
    - $z$
    - $z_p$
    - $y$
*   - Plant Output
    - input to operators
    - performance output
    - output to controller
*   - 
    - $w$
    - $w_p$
    - $u$ 
*   - Plant Input
    - output to operators
    - performance input
    - input to controller -->

The network with state $x^N$ interfacing the operator $F$ and the controller may be described as 
```{math}
\begin{align*}
\mat{c}{x^N_{k+1} \hl z_k \\ z_{p k} \\ y_k} &= \mat{c|cc}{A & B_z & B_{z_p} &  B_u \hl 
C_z & D_{zd} & D_{z w_p} &  D_{zu} \\
C_{z_p} & D_{z_p d} & D_{z_p w_p} &  D_{z_p u} \\
C_y & D_{yd} & D_{y w_p} & D_{yu}} \mat{c}{x_k^N \hl w_k \\ w_{pk} \\ u_k}.
\end{align*}
```



The dimensions of the signals are stored in a cell `n`, such as `n.nz = 4`and `n.nzp = 0`. This partitioning is used to pose the plants. 


## Individual Plants
```{eval-rst}
.. mat:autoclass :: plant.genplant.genplant
    :members:
```


## Cell of Plants (Subsystems)

The `genplant_poly` class is used to store subsystems in the switched systems setting.


```{eval-rst}
.. mat:autoclass :: plant.genplant.genplant_poly
    :members:
```