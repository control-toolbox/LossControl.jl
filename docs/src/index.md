# Loss control regions in optimal control problems
## Introduction

<p>
    The objective of this paper is to address optimal control problems with <i>loss control regions</i>. In that context the state space is partitioned into disjoint sets, referred to as <i>regions</i>, which are classified into two types: <i>control regions</i> and <i>loss control regions</i>. When the state belongs to a control region, the control is permanent (i.e. the control value is authorized to be modified at any time). On the contrary, when the state belongs to a loss control region, the control must remain constant as long as the state belongs to this region. We underline that, in a loss control region, the (constant) value of the control is not predefined: it is a parameter to optimize. Additionally, note that a trajectory can visit several times the same loss control region and, at each visit, the (constant) value of the control can differ.
</p>

![MRI](resources/mri.jpg)



!!! note "2022 survey"

    This background overview is not up-to-date. We refer to for more details.



[^1]: T. Bayen, A. Bouali, L. Bourdin & O. Cots, *Loss control regions in optimal control problems*, Journal of Differential Equations, **12** (2024) 405, 359-397.

## Dependencies

All the numerical simulations to generate this documentation from `MRI.jl` are performed with the following packages.

```@example
using Pkg
Pkg.status()
```