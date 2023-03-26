# Optimal control problems with loss control: applications

Zermelo:

* [permanent case](zermelo/permanent.ipynb) : In this example, we investigate the Zermelo navigation problem under the assumption of permanent control. More precisely, the control input $u(\cdot)$ can be adjusted at any time $t \in [0,8]$.


* [contro loss case 1: one region of control loss](loss1.ipynb) : In this example, we consider the Zermelo navigation problem with a state space divided into two regions: one of type C (control region) and one of type NC (loss of control region). Specifically, the control input $u(\cdot)$ can only be modified when $x(\cdot)$ belongs to the region of type C. However, when $x(\cdot)$ belongs to the region of type NC, the control is constrained to a constant value to be determined.

* [control loss case 1: via augmentaion](zermelo/loss1_aug.ipynb) : In this example, we investigate the same problem as in [contro loss case 1: one region of control loss](loss1.ipynb) only we use the technique of augmentation.

* [contro loss case 2: two regions of control loss](loss2.ipynb) : In this example, we consider the Zermelo navigation problem with a state space divided into three regions: one of type C (control region) and two of type NC (loss of control region). Specifically, the control input $u(\cdot)$ can only be modified when $x(\cdot)$ belongs to the region of type C. However, when $x(\cdot)$ belongs to the regions of type NC, the control is constrained to a constant value to be determined (this constant value may differ between the two regions.)

harmonic oscillator:

* [permanent case](harmonic\ oscillator/permanent.ipynb) : In this example, we investigate the harmonic oscillator problem under the assumption of permanent control. More precisely, the control input $u(\cdot)$ can be adjusted at any time $t$.


* [contro loss case: minimum time problem](loss.ipynb) : In this example, we consider the harmonic oscillator problem, where we minimize the final time, with a state space divided into two regions: one of type C (control region) and one of type NC (loss of control region). Specifically, the control input $u(\cdot)$ can only be modified when $x(\cdot)$ belongs to the region of type C. However, when $x(\cdot)$ belongs to the region of type NC, the control is constrained to a constant value to be determined.


* [contro loss casee 2: minimum energy problem](loss2.ipynb) : In this example, we consider the Zermelo navigation problem, where we minimize the energy, with a state space divided into three regions: one of type C (control region) and two of type NC (loss of control region). Specifically, the control input $u(\cdot)$ can only be modified when $x(\cdot)$ belongs to the region of type C. However, when $x(\cdot)$ belongs to the regions of type NC, the control is constrained to a constant value to be determined (this constant value may differ between the two regions.)
