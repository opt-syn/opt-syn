# Performance Specifications

A performance specification encodes a desired property of the algorithm as a
constraint on the performance channel $(w_p, z_p)$ of the
{doc}`generalized plant <plants/doc_genplant>`. In Analysis the property is verified;
in Synthesis it is the property the returned algorithm is required to meet.

The usage-facing
walkthrough is on the
{doc}`Specifications <../usage/problem_formulation/specs>` page. The channel each
specification acts on is created by the plant's `perf_output_*` and `add_oracle_*` methods
(see {doc}`Generalized Plant <plants/doc_genplant>`).

## Linear Convergence

```{eval-rst}
.. mat:autoclass:: spec.spec_stability
    :members:
```

## Quadratic Performance

```{eval-rst}
.. mat:autoclass:: spec.spec_quad
    :members:
```

## $\ell_2$ Stability (ISS)

```{eval-rst}
.. mat:autoclass:: spec.spec_l2
    :members:
```

## $\ell_2$ Gain

```{eval-rst}
.. mat:autoclass:: spec.spec_e2e
    :members:
```

## Peak-to-Peak

```{eval-rst}
.. mat:autoclass:: spec.spec_p2p
    :members:
```

## Passivity

```{eval-rst}
.. mat:autoclass:: spec.spec_passivity
    :members:
```

## Ergodic Convergence

```{eval-rst}
.. mat:autoclass:: spec.spec_ergodic
    :members:
```