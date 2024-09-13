# A numerical approach for optimal control problems with loss control regions

In optimal control theory, there are several ways for solving numerically an optimal control problem. Direct and indirect methods represent an important class of methods that we will use hereafter. 

Direct methods involve discretizing the state and control variables, simplifying the problem into a finite-dimensional nonlinear optimization problem. On the other hand, indirect methods tackle the problem by solving a boundary value problem, based on the PMP, through the use of a shooting method (see, e.g.[^1],[^2]). 

It is important to note that neither of these methods is fundamentally better than the other. Indeed, each of these methods has its pros and cons. For instance, although the direct method is simple to implement, more robust and less sensitive to the choice of the initial condition, it should be noted that it yields less precise results and can converge to local minima that significantly deviate from the optimal solution. Additionally, this method requires a large amount of memory. On the other hand, the indirect method is known for its extreme precision. However, it is based only on necessary optimality conditions from the PMP and often requires knowledge of the structure of the optimal solution. Moreover, it is quite sensitive to the choice of the initial condition, which must be chosen carefully to ensure convergence.

Often in the literature, one proceeds in two steps. The first step is to implement a direct method to determine the optimal solution's structure and extract the associated adjoint vector. The second step involves constructing an indirect shooting method, where the initial condition is based on the numerical results obtained from the direct method.



## Description of the direct method

For some $\omega_0 \in \mathrm{U}$, some $\varepsilon_0>0$ and $\varepsilon>0$ small enough, we introduce the *regularized problem* given by
```math
\begin{equation}
\begin{array}{lcl}
     \text{minimize}&  & \phi(x(0),x(T)) + \varepsilon_0 \displaystyle\int_0^T v^2(t) \, dt + \displaystyle\int_0^T (1- \Psi_{\varepsilon}(x(t)))\|u(t)-\omega_0\|_{\R^m}^2 \, dt,\\[10pt]
     \text{subject to}& &  (x,\lambda,u,v) \in \mathrm{AC}([0,T],\R^n) \times \mathrm{AC}([0,T],\R^m) \times \mathrm{L}^\infty([0,T],\R^m)\times \mathrm{L}^\infty([0,T],\R), \\[2pt]
     & & \dot{x}(t) = \Psi_{\varepsilon}(x(t)) f(x(t),u(t)) + (1-\Psi_{\varepsilon}(x(t))) f(x(t),\lambda(t)), \quad \text{a.e.\ } t\in [0,T], \\[2pt]
     & & \dot{\lambda}(t) =\Psi_{\varepsilon}(x(t))v(t), \quad \text{a.e.\ } t\in [0,T], \\[2pt]
     &  & g(x(0),x(T)) \in \mathrm{S}, \\[2pt]
     &  & \lambda(t)   \in \mathrm{U}, \quad \text{a.e.\ } t\in [0,T],\\[2pt]
     &  & (u(t),v(t)) \in  \mathrm{U} \times \R , \quad \text{a.e.\ } t\in [0,T],
     \end{array}
\end{equation}
```
where $\Psi_{\epsilon} : \R^n \to \R$ is the regularization of the characteristic function of $\cup_{q_j = 1} \overline{X_j}$ given by 
```math
\Psi_{\varepsilon}(x) := \sum_{q_j = 1} e^{-\frac{1}{2 \varepsilon} \mathrm{d}^2_{j}(x)},
```
for all $x\in \R^n$, where $\mathrm{d}_{j}: \R^n \to \R$ stands for the distance function to the set $\overline{X_j}$ defined by $\mathrm{d}_{j}(x) := \inf_{y \in \overline{X_j}} \|x-y\|_{\R^n}$ for all $x\in \R^n$ and every $j \in \mathcal{J}$. 

## Description of the indirect method

Recall that the direct method has captured the structure of the optimal pair $(x^*,u^*)$. In the indirect method, we address each arc separately. 
We begin by defining the **flow** of the Hamiltonian associated with each arc. To accomplish this, we use the function `Flow` that can be found in the [CTFlows.jl](https://github.com/control-toolbox/CTFlows.jl) package. This latter allows to solve the Hamiltonian system over a given time interval from given initial values of the state and the adjoint vector. This function requires necessary libraries such as  ForwardDiff for calculating gradients and Jacobians and  OrdinaryDiffEq for solving ordinary differential equations. 
    
In the setting of the present paper (including control regions and loss control regions), we distinguish between two types of Hamiltonian flows:
    
- **Hamiltonian flows in control regions** We recall that the Hamiltonian  $H$ associated with the optimal control problem given above is defined by
```math
H(x,u,p) := \langle  p, f(x,u)\rangle_{\R^n},
```
for all $(x,u,p)\in \mathbb{R}^n \times \mathbb{R}^m \times \mathbb{R}^n$. Using the theorem given above and more specifically the Hamiltonian maximization condition, we obtain an expression of the control $u^*$ which can generate a sequence of arcs. Then it remains to define a **pseudo-Hamiltonian** (the pseudo-Hamiltonian stands for the Hamiltonian flow associated with each arc). associated with each arc. Finally we define the flow associated with each arc, which allows the resolution of the boundary value problem on a time interval satisfied by the pair $(x^*,p)$ with an initial condition.
        
- **Hamiltonian flows in loss control regions** Recall that, in loss control regions,  $u^*$ satisfies an averaged Hamiltonian gradient condition (see Theorem above). Here, the difficulty lies in the fact that this condition is given in an integral and implicit form. Therefore, to overcome this difficulty, we first introduce new states  $\lambda$ and  $y$. First, the state $\lambda$ comes from the augmentation technique (we refer to [^3] for more details) to handle the constant value  $u^*_k$, so it satisfies the dynamics  $\dot{\lambda}(t)=0$ and the initial condition  $\lambda(\tau^*_{k-1}) = u^*_k$. Second, the state $y$ satisfies $\dot{y}(t) = 0$ and $y(\tau^*_{k-1})=0$ (and thus $y=0$). Now we define the new Hamiltonian  $\tilde{H}$ as follows:
```math        
\tilde{H}(x,u,y,p) 
            := H(x,u,p) - y \nabla_{u} H(x,u,p) = \langle  p, f(x,u)\rangle_{\R^n} - y \nabla_{u} H(x,u,p),
```
for all  $(x,u,y,p)\in \mathbb{R}^n \times \mathbb{R}^m \times \R \times \mathbb{R}^n $. 
It is important to note that, since $y=0$, we recover the same Hamiltonian  $H$. But, the actual utility of introducing the state  $y$ is that it allows us to rewrite the integral expressed in the averaged Hamiltonian gradient condition as a terminal value of an adjoint vector. This makes it easier to take into account in the shooting function. Here is the justification of this point. First, we define  $p_{y}$ as the solution to the  system
```math
    \left\{\begin{array}{l}
     \dot{p}_{y}(t) = - \nabla_y \tilde{H}(x^*(t),u^*_k,y(t),p(t)), \quad
     \text{for a.e. } t\in [\tau^*_{k-1}, \tau^*_k],  \\
    {p}_{y}(\tau^*_{k-1}) = 0_{\R^n}.
    \end{array}
    \right.
```
Second, since  $u_k^*$ is assumed to be an interior value to  $\mathrm{U}$, we get that
```math
\begin{equation*}
    \int_{\tau^*_{k-1}}^{\tau^*_k} \nabla_{u} H(x^*(t),u^*_k,p(t)) \, dt = {p}_{y}(\tau^*_k)=0.
\end{equation*}
```    
Hence there is no need to compute an integral in order to take into account the averaged Hamiltonian gradient condition. 




















[^1]: S. Aronna, F. Bonnans, P. Martinon, *A shooting algorithm for optimal control problems with singular arcs*, Journal of Optimization Theory and Applications, **158**, pp. 419–459, Springer, 2013.
[^2]: J.F. Bonnans, *The shooting approach to optimal control problems*, IFAC Proceedings Volumes, **46**, no. 11, pp. 281–292, Elsevier, 2013.
[^3]: T. Bayen, A. Bouali, L. Bourdin & O. Cots, *Loss control regions in optimal control problems*, Journal of Differential Equations, **405** (2024), 359-397.
