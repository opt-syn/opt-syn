# Repeated Operator Evaluations

Some optimziation algorithms involve multiple evaluations of operators within a single time step. The Extragradient Method {footcite}`korpelevich1976extragradient`  
```{math}
x_{k+1} = x_{k} - \gamma \nabla f(x_k - \gamma \nabla f(x_k)),
```
can be expressed as the algorithmic interconnection
```{math}
\mat{c}{x_{k+1} \hl z^1_k \\ z^2_k} &= \mat{c|cc}{I & 0 & -\gamma I \hl I & 0 & 0 \\ I & -\gamma I & 0 } \mat{c}{x_k \hl w^1_k \\ w^2_k}, & \mat{c}{w^1_k \\ w^2_k} &= \mat{c}{\nabla f(z^1_k) \\ \nabla f(z^2_k)}.
```

 
The  {attr}`bind` attribute of {class}`opt_system` stores the indexing of repeated operator evaluations.
 For the Extragradient method, this is `bind=[1, 1]`.

To give an example with $s=3$ operators, an algorithm with parameters $\gamma, \lambda > 0$ described by 
```{math}
\mat{c}{x_{k+1} \hl z^1_k \\ z^2_k \\ z^3_k \\ z^4_k} &= 
\mat{c|cccc}{I & -2\gamma \lambda I & -\gamma \lambda I& -2\gamma \lambda I  & -\gamma  \lambda I \hl 
I & -\lambda I & 0 & 0 & 0 \\
I & -3 \lambda I & -\lambda & 0 & 0 \\
I & -\lambda I & 0 & 0 & 0 \\
I & -3 \lambda I & -\lambda I & -2 \lambda I & -\lambda I} \mat{c}{x_k \\ w^1_k \\ w^2_k \\ w^3_k \\ w^4_k}, & \mat{c}{w^1_k \\ w^2_k \\ w^3_k \\ w^4_k} &= \mat{c}{F_1 (z^1_k) \\ F_3 (z^2_k)\\ F_2 (z^3_k) \\ F_3(z^4_k)}
```
evaluates the operator $F_3$ twice in each time step. This algorithm may be modeled using the code
```matlab
Operator_Class = {op1, op2, op3}; %classes for F1, F2, F3
A = [1];
B = [-2, -1, -2, -1] * (-lambda * gamma);
C = [1; 1; 1; 1];
D = [1, 0, 0, 0;
    3, 1, 0, 0;
    1, 0, 0, 0;
    3, 1, 2, 1] * (-lambda);

K = ss(A, B, C, D, 1);

bind = [1, 3, 2, 3]; %ordering of operators in repeated evaluations

sys = opt_system(Operator_Class, [], K);
sys.bind = bind;
```

In the terminology of {footcite}`morin2024frugal`, operator splitting problems with repeated evaluations are 'nonfrugal'.