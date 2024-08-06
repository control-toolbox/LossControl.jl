# Loss control regions in optimal control problems
## Introduction

Often in optimal control problems we consider the case of a permanent control (that can be modified at any instant of time). However, there are situations where the control is lost when the state belongs to a given regions of the state space [^1].

![MRI](mri-resources/mri.jpg)



!!! note "2022 survey"

    This background overview is not up-to-date. We refer to [^18] for more details.



[^1]: T. Bayen, A. Bouali, L. Bourdin & O. Cots, *Loss control regions in optimal control problems*, Journal of Differential Equations, **12** (2024) 405, 359-397.

## Dependencies

All the numerical simulations to generate this documentation from `MRI.jl` are performed with the following packages.

```@example
using Pkg
Pkg.status()
```