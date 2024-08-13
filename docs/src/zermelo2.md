# Zermelo problem, example 1
```math

    \left\{
    \begin{array}{l}
        \displaystyle \min - x_1(8), \\[0.5em]
        \dot{x}_1(t) = x_2(t) + \cos(u(t)), \; \text{for a.e. } t\in [0,8],\\[0.5em]
        \dot{x}_2(t) = \sin(u(t)),  \; \text{for a.e. } t\in [0,8], \\[0.5em]
        u(t) \in [-\frac{\pi}{2}, \frac{\pi}{2}], \; \text{for a.e. } t\in [0,8], \\[0.5em]
        x(0) = 0_{\mathbb{R}^2}, \quad x_2(8) = 4,\\[0.5em]
        \{x \in \mathbb{R}^2 \mid 5 < x_1 < 10 \} \text{ and } \{x \in \mathbb{R}^2 \mid 20 < x_1 < 25 \} \text{ are loss control regions.}
    \end{array}
    \right.
```
# Reformulation for the direct method

```math

    \left\{
    \begin{array}{l}
        \displaystyle \min - x_1(8) + \varepsilon \int_0^8 v^2(t)dt + \int_0^8 f_{NC}(x(t))u^2(t)dt, \\[0.5em]
        \dot{x}_1(t) = f_{C}(x(t)) (x_2(t) + \cos(u(t)))+f_{NC}(x(t)) (x_2(t) + \cos(\lambda(t))), \; \text{for a.e. } t\in [0,8],\\[0.5em]
        \dot{x}_2(t) =  f_{C}(x(t)) \sin(u(t))+f_{NC}(x(t)) \sin(\lambda(t)),  \; \text{for a.e. } t\in [0,8], \\[0.5em]
        \dot{\lambda}(t) = f_{C}(x(t))v(t),  \; \text{for a.e. } t\in [0,8], \\[0.5em]
        u(t)\in [-\frac{\pi}{2}, \frac{\pi}{2}] , \; \text{for a.e. } t\in [0,8], \\[0.5em]
        x(0) = 0_{\mathbb{R}^2}, \quad x_2(8) = 4,\\[0.5em]
        \{x \in \mathbb{R}^2 \mid 5 < x_1 < 10 \} \text{ and } \{x \in \mathbb{R}^2 \mid 20 < x_1 < 25 \} \text{ are loss control regions.}
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
# Hyperbolic tangent function for smooth indicator
function smooth_indicator_tanh(x, a, b, ε)
    return 0.5 * (tanh((x - a) / ε) - tanh((x - b) / ε))
end

# Parameters for the intervals [5, 10] and [20, 25]
a1, b1 = 5.0, 10.0
a2, b2 = 20.0, 25.0
ε2   = 1e-3
tf  = 8

# Regularization parameter
ε1 = 0.05 # You can adjust this to control the smoothness

# Define x range
x = range(0, 30, length=1000)

# Control regions/loss control regions indicator functions
fNC(x) = smooth_indicator_tanh.(x, a1, b1, ε1) .+ smooth_indicator_tanh.(x, a2, b2, ε1)
fC(x)  = 1 .- fNC(x)

# Plots
plot(x, fNC, label="fNC")
plot!(x, fC, label="fC")
xlabel!("x")
ylabel!("y")
title!("Smooth Indicator Function for [5,10] U [20,25]")

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
     q̇(t) == [fNC(x1(t))*(x2(t) + cos(λ(t))) + fC(x1(t))*(x2(t) + cos(u(t))),
             fNC(x1(t))*sin(λ(t)) + fC(x1(t))*sin(u(t)),
             fC(x1(t))*v(t),
             (v(t))^2,
             fNC(x1(t))*(u(t))^2]

    #cost function        
    -x1(tf) + ε2*xv(tf) + xu(tf) → min    
end;
```


```@example main
N = 400
sol = solve(ocp; grid_size=N);
```


```@example main
plot(sol; layout=:group, size=(800, 300))
```
