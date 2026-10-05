# Performance specifications

A performance specification `spec` encodes a desired property of the algorithm as a
constraint on the performance channel $(w_p, z_p)$ of the
{doc}`generalized plant <plants/doc_genplant>`. In analysis, `spec` represents the  property to be  verified.
In synthesis, `spec` is the property that the returned algorithm should meet.

The usage-facing
walkthrough is on the
{doc}`specifications <../usage/problem_formulation/specs>` page. The channel each
specification acts on is created by the plant's `perf_output_*` and `add_oracle_*` methods
(see {doc}`generalized plant <plants/doc_genplant>`).

## Linear convergence

```{eval-rst}
.. mat:autoclass:: spec.spec_stability
    :members:
```

## Quadratic performance

```{eval-rst}
.. mat:autoclass:: spec.spec_quad
    :members:
```

## $\ell_2$ stability (ISS)

```{eval-rst}
.. mat:autoclass:: spec.spec_l2
    :members:
```

## $\ell_2$ gain

```{eval-rst}
.. mat:autoclass:: spec.spec_e2e
    :members:
```

<!-- ## Peak-to-Peak

```{eval-rst}
.. mat:autoclass:: spec.spec_p2p
    :members:
``` -->

## Passivity

```{eval-rst}
.. mat:autoclass:: spec.spec_passivity
    :members:
```

## Stochastic sensitivity

```{eval-rst}
.. mat:autoclass:: spec.spec_h2
    :members:
```

<!-- ## Ergodic Convergence

```{eval-rst}
.. mat:autoclass:: spec.spec_ergodic
    :members:
``` -->