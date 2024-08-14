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
