# Direct method
```math
    \left\{
    \begin{array}{l}
        \displaystyle \min - x_1(8), \\[0.5em]
        \dot{x}_1(t) = x_2(t) + \cos(u(t)), \; \text{for a.e. } t\in [0,8],\\[0.5em]
        \dot{x}_2(t) = \sin(u(t)),  \; \text{for a.e. } t\in [0,8], \\[0.5em]
        u(t) \in [-\frac{\pi}{2}, \frac{\pi}{2}], \; \text{for a.e. } t\in [0,8], \\[0.5em]
        x(0) = 0_{\mathbb{R}^2}, \quad x_2(8) = 4,\\[0.5em]
        \{x \in \mathbb{R}^2 \mid 0.5 < x_2 < 3.5 \} \text{ is a loss control region.}
    \end{array}
    \right.
```

```math
    \left\{
    \begin{array}{l}
        \displaystyle \min - x_1(8) + \epsilon\int_0^8 v^2(t)dt + \int_0^1 f_{NC}(x(t))u^2(t)dt, \\[0.5em]
        \dot{x}_1(t) = f_{C}(x(t))(x_2(t) + \cos(u(t))) + f_{NC}(x_2(t) + \cos(\lambda(t))), \; \text{for a.e. } t\in [0,8],\\[0.5em]
        \dot{x}_2(t) = f_{C}(x(t))\sin(u(t)) + f_{NC}(x(t))\sin(\lambda(t)),  \; \text{for a.e. } t\in [0,8], \\[0.5em]
        \dot{\lambda}(t) = f_{C}(x(t))v^2(t),  \; \text{for a.e. } t\in [0,8], \\[0.5em]

        u(t) \in [-\frac{\pi}{2}, \frac{\pi}{2}], \; \text{for a.e. } t\in [0,8], \\[0.5em]
        x(0) = 0_{\mathbb{R}^2}, \quad x_2(8) = 4,\\[0.5em]
        \{x \in \mathbb{R}^2 \mid 0.5 < x_2 < 3.5 \} \text{ is a loss control region.}
    \end{array}
    \right.
```


```@example main
using JuMP  
using Ipopt
using Plots
using Plots.PlotMeasures
using LaTeXStrings
```


```@example main
using OptimalControl
using NLPModelsIpopt
```


```@example main
#hyperbolic tangent function
function smooth_indicator_tanh(x, a, b, ε)
    return 0.5 * (tanh((x - a) / ε) - tanh((x - b) / ε))
end

#parameters
a, b = 0.5, 3.5  
ε = 0.01  
x = range(0, 4, length=1000)
ε   = 1e-3
tf  = 8

#control regions/loss control regions indidcator functions
fNC(x)        = smooth_indicator_tanh.(x, a, b, ε) 
fC(x)         = 1-smooth_indicator_tanh.(x, a, b, ε)

#plots
plot(x, fNC, label="fNC")
plot!(x, fC, label="fC")
xlabel!("x")
ylabel!("y")
title!("Smooth Indicator Function Approximations");
```


```@example main
@def ocp begin

    t ∈ [ 0., tf ],                  time

    q = [ x1, x2, λ, xu, xv ] ∈ R^5, state

    ω = [u, v] ∈ R^2,                control

    #initial conditions
    x1(0) == 0
    x2(0) == 0
    xu(0) == 0
    xv(0) == 0

    #final condition
    x2(tf) == 4

    #control constraint
    -π/2  ≤  u(t)  ≤ π/2

    #state constraints
    -π/2  ≤  λ(t)  ≤ π/2

    #hybrid control system
     q̇(t) == [fNC(x2(t))*(x2(t) + cos(λ(t))) + fC(x2(t))*(x2(t) + cos(u(t))),
             fNC(x2(t))*sin(λ(t)) + fC(x2(t))*sin(u(t)),
             fC(x2(t))*v(t),
             (v(t))^2,
             fNC(x2(t))*(u(t))^2]

    #cost function        
    -x1(tf) + ε*xv(tf) + xu(tf) → min    
end;
```


```@example main
N = 500
sol = solve(ocp; grid_size=N);
```




```@example main
plot(sol; layout=:group, size=(800, 300))
```
