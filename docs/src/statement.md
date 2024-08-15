# Statement of the problem

Let $n$, $m$, $\ell \in \mathbb{N}^*$ be three positive integers and $T > 0$ be a positive real number. In this section, we consider a partition of the state space given by

```math
\mathbb{R}^n = \bigcup_{j \in \mathcal{J}} \overline{X_j},
```

where $\mathcal{J}$ is a (possibly infinite) family of indexes and the nonempty open subsets $X_j \subset \mathbb{R}^n$, called *regions*, are disjoint. We introduce an indexation $q_j \in \{0, 1\}$ allowing us to separate *control regions* and *loss control regions* (see Introduction for details) as follows

```math
q_j := \begin{cases}
1 \text{ if } X_j \text{ is a control region,} \\
0 \text{ if } X_j \text{ is a loss control region,}
\end{cases}
```

for all $j \in \mathcal{J}$.

Our aim in this section is to derive first-order necessary optimality conditions in a Pontryagin form for the optimal control problem with *loss control regions* given by

```math
\begin{aligned}
\text{minimize} & \quad \phi(x(0),x(T)), \\
\text{subject to} & \quad (x,u) \in \mathrm{AC}([0,T],\mathbb{R}^n) \times \mathrm{L}^\infty([0,T],\mathbb{R}^m), \\
& \quad \dot{x}(t) = f(x(t), u(t)), \quad \text{a.e.\ } t \in [0,T], \\
& \quad g(x(0),x(T)) \in \mathrm{S}, \\
& \quad u(t) \in \mathrm{U}, \quad \text{a.e.\ } t \in [0,T], \\
& \quad u \text{ is constant when $x$ is in a loss control region},
\end{aligned}
```

where the Mayer cost function $\phi: \mathbb{R}^n \times \mathbb{R}^n \to \mathbb{R}$, the dynamics $f: \mathbb{R}^n \times \mathbb{R}^m \to \mathbb{R}^n$, and the constraint function $g: \mathbb{R}^n \times \mathbb{R}^n \to \mathbb{R}^\ell$ are of class $\mathrm{C}^1$, and where both subsets $\mathrm{S} \subset \mathbb{R}^\ell$ and $\mathrm{U} \subset \mathbb{R}^m$ are nonempty closed convex subsets.

## Regular solution to the control system with loss control regions

First, let us provide the definition of a solution to the following control system

```math
\begin{equation*}
\left\lbrace 
\begin{array}{l}
\dot{x}(t) = f(x(t),u(t)),  \quad \text{for a.e. } t \in [0,T],  \\[2pt]
u \text{ is constant when } x \text{ is in a loss control region}.
\end{array}
\right.
\quad \text{(CS)}
\end{equation*}
```

!!! note "Definition (Solution to (CS))"
    A pair $(x,u) \in \mathrm{AC}([0,T],\R^n) \times \mathrm{L}^\infty([0,T],\R^m)$ is said to be a solution to (CS) if there exist a finite number $\mathbb{N}^*$ and a partition $\mathbb{T} = \{\tau_k\}_{k=0,\ldots,N}$ of the interval $[0,T]$ such that:               
    - It holds that
    ```math
        \forall k \in \{ 1,\ldots,N \}, \quad \exists j(k) \in \mathcal{J}, \quad \forall t \in (\tau_{k-1},\tau_k), \quad x(t) \in X_{j(k)}, 
    ```
    where $j(k) \neq j(k-1)$ for all $k \in \{ 2,\ldots,N \}$. The sequence $\{j(1),\ldots, j(N)\}$ is called the *switching sequence*.
    - It holds that $x(0) \in X_{j(1)}$ and $x(T) \in X_{j(N)}$.
    - For all $k \in \{1, \ldots, N\}$ such that $q_{j(k)}=0$, the control $u$ is constant over $(\tau_{k-1}, \tau_k)$ (the constant value being denoted by $u_k$ in the sequel).
    - It holds that $\dot x(t)=f(x(t),u(t))$ for almost every $t \in [0,T]$.
    The times $\tau_k$ for $k \in \{1,\ldots,N-1\}$, called  \textit{crossing times}, correspond to the instants at which the trajectory $x$ goes from the region $X_{j(k)}$ to the region $X_{j(k+1)}$, and thus $x(\tau_k) \in \partial X_{j(k)} \cap \partial X_{j(k+1)}$. 

The PMP with loss control regions is based on some regularity assumptions made on the optimal pair of the optimal control problem with loss control regions given above at each of its crossing times. These hypotheses are made more precise in the next definition.

!!! note "Definition (Regular solution to (CS))"
    Following the notations introduced in Definition given above, a solution $(x,u) \in \mathrm{AC}([0,T],\R^n) \times \mathrm{L}^\infty([0,T],\R^m)$ to(CS), associated with a finite number $N \in \mathbb{N}^*$, a partition $\mathbb{T} = \{ \tau_k \}_{k=0,\ldots,N}$ and a switching sequence $\{j(1), \ldots, j(N)\}$, is said to be *regular* if the following conditions are both satisfied:
    - At each crossing time $\tau_k$, there exists a $\mathrm{C}^1$ function $F_k : \R^n \to \R$ such that 
    ```math
        \begin{equation}
            \exists \nu_k > 0, \quad
            \forall z \in  \overline{\mathrm{B}}_{\R^n}(x(\tau_k),\nu_k), \quad  
            \left\{
            \begin{array}{rcl} 
                z \in X_{j(k)} & \Leftrightarrow & F_k(z)<0, \\[2pt]
                z \in \partial X_{j(k)}\cap \partial X_{j(k+1)} & \Leftrightarrow & F_k(z)=0, \\[2pt]
                z \in X_{j(k+1)} & \Leftrightarrow & F_k(z)>0.
            \end{array}
            \right.
        \end{equation}
    ```
    In particular it holds that $F_k(x(\tau_k))=0$.
    - At each crossing time $\tau_k$, there exists $\alpha_k>0$ and $\beta_k > 0$ such that the *transverse condition*
    ```math
    \begin{equation}
    \langle  \nabla F_k(x(\tau_k)) , f(x(\tau_k),u(t)) \rangle_{\R^n} \geq \beta_k,  \quad \text{a.e.\ } t \in [\tau_k- \alpha_k, \tau_k+\alpha_k],
    \end{equation}
    ```
    is satisfied.

## Pontryagin maximum principle with loss control regions

The Hamiltonian $H : \R^n \times \R^m \times \R^n \to \R$ associated with Problem \eqref{P} is defined by
```math
 H(x,u,p) := \langle  p , f(x,u) \rangle_{\R^n}
``` 
for all $(x,u,p) \in \R^n \times \R^m \times \R^n$. We are now in a position to state the main result of this section.

!!! tip "Theorem"
    If $(x^*,u^*)\in \mathrm{AC}([0,T],\R^n) \times \mathrm{L}^{\infty}([0,T],\R^m)$ is a global solution to the above problem, that is moreover a regular solution to (CS), associated with a finite number $N \in \mathbb{N}^*$, a partition $\mathbb{T}^* = \{ \tau^*_k \}_{k=0,\ldots,N}$ and a switching sequence $\{j(1),\ldots,j(N)\}$, and such that $g$ is submersive at $(x^*(0),x^*(T))$, then there exists a nontrivial pair $(p,p^0) \in \mathrm{PAC}_{\mathbb{T}^*}([0,T],\R^n) \times \R_+$ satisfying:
    -  the *Hamiltonian system*
    ```math
    \begin{equation*}
    \dot{x^*}(t) = \nabla_p H(x^*(t),u^*(t),p(t)) \quad \text{and} \quad
    -\dot{p}(t) = \nabla_x H(x^*(t),u^*(t),p(t)),
    \end{equation*}
    ```
    for almost every $t \in [0,T]$;
    - the *transversality condition*
    ```math
    \begin{equation*}
    \left( \begin{array}{c}
    p(0) \\[5pt]
    -p(T) 
    \end{array} \right)
    = p^0 \nabla \phi (x^*(0),x^*(T)) + \nabla g (x^*(0),x^*(T)) \xi,
    \end{equation*}
    ```
    for some $\xi \in \mathrm{N}_\mathrm{S}[ g (x^*(0),x^*(T)) ]$;
    - the *discontinuity condition*
    ```math
    \begin{equation*}
    p^+(\tau^*_k) - p^-(\tau^*_k) =  \sigma_k
    \nabla F^*_k(x^*(\tau^*_k)),
    \end{equation*}
    ```
    for some $\sigma_k \in \R$ and for all $k \in \{ 1,\ldots , N-1\}$; 
    - the *Hamiltonian maximization condition*
    ```math
    \begin{equation*}
        u^*(t) \in \mathrm{argmax}_{\omega \in \mathrm{U}}H(x^*(t),\omega,p(t)),
    \end{equation*}
    ```
    for almost every $t\in (\tau^*_{k-1},\tau^*_k)$ and all $k \in \{1,\ldots,N\}$ such that $q_{j(k)}=1$;
    - the *averaged Hamiltonian gradient condition*
    ```math
    \begin{equation*}
       \int_{\tau^*_{k-1}}^{\tau^*_k} \nabla_{u} H(x^*(t),u^*_k,p(t))\, dt \in \mathrm{N}_{\mathrm{U}}[u^*_k],
    \end{equation*}
    ```
    for all $k \in \{1,\ldots,N\}$ such that $q_{j(k)}=0$;
    - the *Hamiltonian constancy condition* 
    ```math
    H(x^*(t),u^*(t),p(t)) = c,
    ```
    for almost every $t\in [0,T]$, for some $c \in \R$.

