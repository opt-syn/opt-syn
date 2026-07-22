# Periodic Systems


A Periodic system has a representation
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k } \mat{c}{x_k \\ w_k}.
```

in which there exists a period $h$ such that
```{math}
 \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k }  = \mat{c|c}{\Acl_{k+h} & \Bcl_{k+h} \hl \Ccl_{k+h} & \Dcl_{k+h} } 
```

These methods are relevant when the network and/or controller are periodic.