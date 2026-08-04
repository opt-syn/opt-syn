# Set-Valued Maps

We support set-valued maps $w \in F(z)$ that satisfy properties
- maximal monotonicity,
- $\mu$-strong-monotonicity $\mu > 0$,
- $\mu$-hypo-monotonicity with $\mu > 0$,
- $\beta$-cocoercivity with $\beta>0$,
- $L$-Lipschitzness with $L > 0$,
- $R$-Inverse-Lipschitzness with $R > 0$. 


The properties are stored in a `prop` cell.



 As an example, an operator that is $\mu$-strongly monotone and $\beta$-cocoercive can be declared using the property structure `prop = {'monotone', mu, 'cocoercive', beta}`.  

The `order` supplied to Analysis is a single integer (number of lags).

```{eval-rst}
.. mat:autoclass :: operator.op_gen  
    :members:
```