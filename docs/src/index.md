# Loss control regions in optimal control problems
## Introduction

<p>
    The objective of this paper is to address optimal control problems with <i>loss control regions</i>. In that context the state space is partitioned into disjoint sets, referred to as <i>regions</i>, which are classified into two types: <i>control regions</i> and <i>loss control regions</i>. When the state belongs to a control region, the control is permanent (i.e. the control value is authorized to be modified at any time). On the contrary, when the state belongs to a loss control region, the control must remain constant as long as the state belongs to this region. We underline that, in a loss control region, the (constant) value of the control is not predefined: it is a parameter to optimize. Additionally, note that a trajectory can visit several times the same loss control region and, at each visit, the (constant) value of the control can differ.
</p>

<p>
    To address optimal control problems with loss control regions, we pursue the approach initiated in our previous works <a href="#bayen2022hybrid">[1]</a>, <a href="#proceeding">[2]</a> which consists in duplicating the control: one permanent control that intervenes only in control regions, and one <i>regionally switching parameter</i><sup>1</sup> that intervenes only in loss control regions. With this approach, we are confronted with a discontinuous dynamics, and thus our framework falls into the domain of <i>hybrid</i> optimal control theory.
</p>

<p>
    In that field, the so-called <i>hybrid maximum principle</i> (in short, HMP) extends the PMP to various hybrid settings (i.e. to discontinuities of various natures). In the HMP, the adjoint vector is usually piecewise absolutely continuous, admitting a discontinuity jump at each time the dynamics changes discontinuously. We refer to <a href="#caines2006">[3]</a>, <a href="#DMITRUK2008">[4]</a>, <a href="#piccoli2005">[5]</a>, <a href="#pakniyat2020hybrid">[6]</a>, <a href="#caines2003">[7]</a>, <a href="#shaikh2007hybrid">[8]</a>, <a href="#sussmann1999maximum">[9]</a> and references therein.
</p>

<p>
    To be specific, our approach leads to a dynamics that changes discontinuously only according to the state position in the state space: we speak of a <i>spatially heterogeneous dynamics</i> and this setting corresponds to the spirit of previous works such as <a href="#trelat2016">[10]</a>, <a href="#trelat2011">[11]</a> (in which transversal crossing assumptions are made to handle the boundary crossings of the optimal trajectory). In contrast, one can find in the literature several contexts where the change of dynamics depends on the time variable only, i.e. one fixes in advance a certain number of instants (fixed or free) at which the dynamics changes. These discontinuities may be controlled or not (see, e.g., <a href="#caines2006">[3]</a>, <a href="#DMITRUK2008">[4]</a>). In such a case, we would rather speak of <i>temporal discontinuities</i>. We emphasize that this type of discontinuities does not appear all along this paper in which only <i>spatial discontinuities</i> are involved.
</p>

<p>
    <small>1. A <i>regionally switching parameter</i> is a parameter that can change its value when the trajectory moves from one region to another. To our best knowledge, this concept has never been considered in the literature until our previous work <a href="#bayen2022hybrid">[1]</a>.</small>
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