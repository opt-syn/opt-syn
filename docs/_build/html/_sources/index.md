---

tocdepth: 1

---

# {{osyn}}



Analysis and Synthesis of First-Order Algorithms (MATLAB)

---

{{osyn}}: <span style="color:#5D93BF">optimization</span> (algorithm) synthesis.

## Overview


{{osyn}}  analyzes and synthesizes first-order algorithms by using methods from robust control theory. The target performance criteria for Analysis includes worst-case linear convergence rates, and gain bounds for the amplification of errors. The Synthesis process tries to design an optimization algorithm satisfying these desired constraints. 

The Analysis and Synthesis problems are posed as convex problems with  Linear Matrix Inequality constraints, and are solved using [LMILab](https://www.mathworks.com/help/robust/ug/introduction.html).  

These  
Analysis and Synthesis tasks may be performed for algorithms arising in dynamic environments,  including cases  with constant or time-varying delays, channel memory, and cross-talk.



## Get Started

For installation and  examples of algorithm  analysis and synthesis workflows, see  {doc}`Get Started <pages/get_started>`.

## Contributors

- [Jared Miller](https://jarmill.github.io/): Creator and maintainer
- [Fabian Jakob](https://www.ist.uni-stuttgart.de/institute/team/Jakob-00004/): Creator
- Manuel Zobel: Creator
- [Carsten Scherer](https://www.imng.uni-stuttgart.de/institute/team/Scherer-00006/): Creator
- [Andrea Iannelli](https://www.ist.uni-stuttgart.de/institute/team/Iannelli/): Creator


Contact [Jared Miller](https://jarmill.github.io/) for questions, issues, comments, and suggestions.

```{toctree}
:maxdepth: 2
:hidden:

Home <self>
Get Started <pages/get_started>
How it Works <pages/how_it_works/index_how_it_works>
Usage <pages/usage/index_usage>
Examples <pages/examples>
Documentation <pages/documentation/index_documentation>
Changelog <pages/changelog/changelog>
Resources <pages/resources>
```


