


## Harmonic oscillator problem

```math
    \left\{
    \begin{array}{l}
        \displaystyle \min  T, \\[0.5em]
        \dot{x}_1(t) = x_2(t), \; t\in [0,T]\\[0.5em]
        \dot{x}_2(t) = u(t)-x_1(t), t\in [0,T]  \\[0.5em]
        u(t) \in [-1, 1], \; t\in [0,T]\\[0.5em]
        x(0) = (4.2,0) , \quad x(T) = 0_{\mathrm{R}^2}, \\[0.5em]
        \{x \mid x_2 < 0\} \text{ is a control loss reigon.}
    \end{array}
    \right.
```

## Reformulation for the direct method

```math
    \left\{
    \begin{array}{l}
        \displaystyle \min  T + \varepsilon \int_0^T v^2(t)dt + \int_0^T f_{NC}(x_2(t))u^2(t)dt, \\[0.5em]
        \dot{x}_1(t) = x_2(t), \; t\in [0,T]\\[0.5em]
        \dot{x}_2(t) =f_{C}(x_2(t))(u(t) - x_1(t)) + f_{NC}(x_2(t))(\lambda(t) - x_1(t))
        , t\in [0,T]  \\[0.5em]
        \dot{\lambda}(t) = f_C(x_2(t))v(t),\; t\in [0,T]\\[0.5em]
        u(t) \in [-1, 1] \; t\in [0,T]\\[0.5em]
        x(0) = (4.2,0) , \quad x(T) = 0_{\mathrm{R}^2}.
    \end{array}
    \right.
```

```@example main
    using JuMP  
    using Ipopt
    using Plots
    using Plots.PlotMeasures
    using LaTeXStrings
    using OptimalControl
    using NLPModelsIpopt
    include("smooth.jl");
```


```@example main
a  = 0.0 
ε1 = 0.02
fNC(x) = fNC_unboundedminus(x,a,ε1)
plot(fNC,-1., 1, label="fNC")
```


```@example main
@def ocp begin
        
    ε = 1e-3

    tf ∈ R,                          variable

    t ∈ [ 0., tf ],                  time

    q = [ x1, x2, λ, xu, xv ] ∈ R^5, state

    ω = [u, v] ∈ R^2,                control

    tf ≥ 0.
    
    #initial conditions
    x1(0) == 2.5
    x2(0) == 4
    xu(0) == 0
    xv(0) == 0

    #final condition
    x1(tf) ==  1e-6
    x2(tf) == -1e-6

    #control constraint
    -1. ≤  u(t)  ≤ 1.

    #state constraint
    -1 ≤  λ(t) ≤ 1,             (1)
    -5 ≤ x1(t) ≤ 5,             (2)
    -5 ≤ x2(t) ≤ 5,             (3)

    #hybrid control system
    q̇(t) == [x2(t), (1-fNC(x2(t)))*u(t) + fNC(x2(t))*λ(t) - x1(t), (1-fNC(x2(t)))*v(t), (v(t))^2, fNC(x2(t))*(u(t))^2]

    #cost function        
    tf + ε*xv(tf) + xu(tf) → min    
end
```


```@example main
    sol = solve(ocp; init = (state = t -> [0.1, 0.1, 1, 0, 0], control =[-1, 0], variable =15), max_iter=650)
```


```@example main
    plot(sol; layout=:group, size=(800, 300))
```