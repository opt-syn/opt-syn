# Optimization Algorithm Synthesis

This is a collection of m-files for the design of optimization algorithms and extremum controllers.

The related theory is described in the paper

C.W. Scherer, Ch. Ebenbauer, T. Holicki, **Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial**, 62nd IEEE Conference on Decision and Control.

[An extended version of this paper is also available on arXiv.](https://doi.org/10.48550/arXiv.2306.00565)


Slight modifications to this work were done by J. Miller. 

All intellectual property remains with the original authors: C.W. Scherer, Ch. Ebenbauer, T. Holicki

## List of files to be available in the current Matlab path

**Files to generate results in paper**

Demo_rates: Generates Fig.  7\
Demo_algos: Generates Figs. 8-10


**Main files for synthesis**

syzf.m\
syzfb.m\
syzfco.m

**Auxiliary file for simulation of algorithm**

algsim.m

**General auxiliary files**

jordanb.m\
rhotrafo.m\
sym2fun.m

## Toolboxes required to run the code

Control System Toolbox\
Robust Control Toolbox\
Yalmip

## Acknowledgements of funding 

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy - EXC 2075 – 390740016. We acknowledge the support by the Stuttgart Center for Simulation Science (SimTech).
