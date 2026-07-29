---

tocdepth: 1

---

# {{osyn}}



Analysis and Synthesis of First-Order Algorithms (MATLAB)

---

{{osyn}}: <span style="color:#5D93BF">optimization</span> (algorithm) synthesis.

## Overview


{{osyn}} uses methods analyzes and synthesizes first-order optimization algorithms from robust control theory. The target performance criteria for Analysis includes worst-case linear convergence rates and gain bounds for the amplification of errors. The Synthesis process tries to design an optimization algorithm satisfying these desired constraints. 
The Analysis and Synthesis problems are posed as convex problems with  Linear Matrix Inequality constraints, and are solved using [LMILab](https://www.mathworks.com/help/robust/ug/introduction.html).  These algorithms 
Analysis and Synthesis may be performed for algorithms arising in networked environments,  including cases with channel memory, constant or time-varying delays, and cross-talk.



## Get Started

For installation and  examples of algorithm  analysis and synthesis workflows, see  {doc}`Get Started <pages/get_started>`.

## Contributors

- [Jared Miller](https://jarmill.github.io/): Creator and maintainer.
- [Fabian Jakob](https://www.ist.uni-stuttgart.de/institute/team/Jakob-00004/): Creator.
- [Carsten Scherer](https://www.imng.uni-stuttgart.de/institute/team/Scherer-00006/): Creator.
- [Andrea Iannelli](https://www.ist.uni-stuttgart.de/institute/team/Iannelli/): Creator.


```{toctree}
:maxdepth: 2
:hidden:

Home <self>
Get Started <pages/get_started>
How it Works <pages/how_it_works/index_how_it_works>
Usage <pages/usage/index_usage>
Examples <pages/examples>
API Reference <pages/documentation/index_documentation>
Changelog <pages/changelog/changelog>
Resources <pages/resources>
```


