# Loss control regions in optimal control problems
## Introduction

**General context.** *Optimal control theory* studies controlled systems to achieve desired targets with minimal cost. The *Pontryagin maximum principle* (PMP) provides necessary conditions for optimality, ensuring an *adjoint vector* (or *costate*) meets the *Hamiltonian maximization condition*.

Typically, optimal control involves *permanent control*, allowing modification of the control function at each time instant. However, practical constraints can lead to *nonpermanent control*. For instance, digital controls result in *sampled-data control* with discrete changes. In aerospace, *eclipse constraints* limit control for solar-powered satellites in a shadow region where the control is reduced to zero. Hence, it is desirable to keep the system outside these regions.

 ![MRI](resources/mri.jpg) 

**Objective and approach.** Here, we address optimal control problems with *loss control regions*[^1], where the state space is divided into *control regions* and *loss control regions*. In control regions, control can change at any time, while in loss control regions, control must remain constant, though its value is to be optimized and can vary with each visit.

We extend our previous work by using a permanent control for control regions and a *regionally switching parameter* for loss control regions. This leads to a discontinuous dynamics framework, fitting into *spatially hybrid optimal control theory*. The *hybrid maximum principle* (HMP) extends the PMP to hybrid settings, with a piecewise absolutely continuous adjoint vector.

**Numerical contribution.** In this note we illustrate a two-step numerical method for optimal control problems with loss control regions. First, a direct numerical approach is applied to a *regularized* problem to manage discontinuities and outline the optimal trajectory's structure. Second, this helps initialize an indirect numerical method for the original problem, using the PMP from Theorem~1. The method incorporates the averaged Hamiltonian gradient condition and adjoint vector discontinuities to define an appropriate shooting function, adding to classical terms for non-hybrid optimal control problems (see ref).

!!! note "Contents"
    - Provide a statement of the PMP with loss control regions.
    - Provide a direct method for solving optimal control problems with loss control regions (based on a regularization technique).
    - Provide an indirect method (shooting method) for solving optimal control problems with loss control regions using the PMP with loss control regions.




[^1]: T. Bayen, A. Bouali, L. Bourdin & O. Cots, *Loss control regions in optimal control problems*, Journal of Differential Equations, **12** (2024) 405, 359-397.

## Dependencies

All the numerical simulations to generate this documentation from `MRI.jl` are performed with the following packages.

```@example
using Pkg
Pkg.status()
```