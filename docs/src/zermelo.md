# Zermelo problem, example 1
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
## Reformulation for the direct method

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
using OptimalControl
using NLPModelsIpopt
include("smooth.jl");
nothing # hide
```




```@example main
I = [(0.5, 3.5)]
ε1 = 0.01  
fNC(x) = fNC_bounded(x,I,ε1)
plot(fNC,0., 5, label="fNC")
```


```@example main
@def ocp begin

    ε   = 1e-3
    tf  = 8


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
     q̇(t) == [fNC(x2(t))*(x2(t) + cos(λ(t))) + (1-fNC(x2(t)))*(x2(t) + cos(u(t))),
             fNC(x2(t))*sin(λ(t)) +(1-fNC(x2(t)))*sin(u(t)),
             (1-fNC(x2(t)))*v(t),
             (v(t))^2,
             fNC(x2(t))*(u(t))^2]

    #cost function        
    -x1(tf) + ε*xv(tf) + xu(tf) → min    
end;
nothing # hide
```


```@example main
N = 500
sol = solve(ocp; grid_size=N); 
nothing # hide
```



```@example main
plot(sol; layout=:group, size=(800, 300))
```

```@example main
tt    = sol.times
tf    = 8
x1(t) = sol.state(t)[1]
x2(t) = sol.state(t)[2]
λ(t)  = sol.state(t)[3]
u(t) = sol.control(t)[1]
nothing # hide
```

```@example main
plot(x1, x2, 0, tf, label="optimal trajectory", color="blue", linewidth=2)
plot!([0, 31], [0.5, 0.5], color=:black, label = false, linewidth=2)
plot!([0, 31], [3.5, 3.5], color=:black, label = false, linewidth=2)
```


```@example main
plot(tt, u, label="optimal control", color="red", linewidth=2)
plot!(tt, λ, label="state λ", color="green", linewidth=2)
```

```@example main
# Find the first crossing time
t1_index = findfirst(t -> x2(t) ≥ 0.5, tt)

# If t1 is found, find the next crossing time
if t1_index !== nothing
    t1       = tt[t1_index]
    t2_index = findfirst(t ->  x2(t) ≥ 3.5, tt[t1_index+1:end])
    t2_index = t2_index !== nothing ? t2_index + t1_index : nothing
    t2 = t2_index !== nothing ? tt[t2_index] : "No such t2 found"
else
    t1 = "No such t1 found"
    t2 = "No such t2 found"
end
println("first crossing time: ",  t1)
println("second crossing time: ", t2)
```

# Zermelo problem, example 2

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
## Reformulation for the direct method

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
I = [(5, 10), (20, 25)]
ε1 = 0.05  
fNC(x) = fNC_bounded(x,I,ε1)
plot(fNC,0., 30, label="fNC")
```


```@example main
@def ocp1 begin

    ϵ  = 1e-3
    
    tf = 8

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
     q̇(t) == [fNC(x1(t))*(x2(t) + cos(λ(t))) + (1-fNC(x1(t)))*(x2(t) + cos(u(t))),
             fNC(x1(t))*sin(λ(t)) + (1-fNC(x1(t)))*sin(u(t)),
             (1-fNC(x1(t)))*v(t),
             (v(t))^2,
             fNC(x1(t))*(u(t))^2]

    #cost function        
    -x1(tf) + ϵ*xv(tf) + xu(tf) → min    
end;
nothing # hide
```


```@example main
N = 400
sol1 = solve(ocp1; grid_size=N);
nothing # hide
```


```@example main
plot(sol; layout=:group, size=(800, 300))
```

```@example main
tt    = sol1.times
tf    = 8
x1(t) = sol1.state(t)[1]
x2(t) = sol1.state(t)[2]
λ(t)  = sol1.state(t)[3]
u(t)  = sol1.control(t)[1];
nothing # hide
```

```@example main 
# Plot the optimal trajectory
plot(x1, x2, 0, tf, label="optimal trajectory", color="blue", linewidth=2)
plot!([5, 5], [0, 5], color=:black, label=false, linewidth=2)
plot!([10, 10],  [0, 5], color=:black, label=false, linewidth=2)
plot!([20, 20],  [0, 5], color=:black, label=false, linewidth=2)
plot!([25, 25],  [0, 5], color=:black, label=false, linewidth=2)
```

```@example main
plot(tt, u, label="optimal control", color="red", linewidth=2)
plot!(tt, λ, label="state λ", color="green", linewidth=2)
```

```@example main
# Find the crossing times based on conditions for x1
t1_index = findfirst(t -> x1(t) > 5, tt)
t2_index = nothing
t3_index = nothing
t4_index = nothing

# If t1 is found, find the next crossing times
if t1_index !== nothing
    t2_index = findfirst(t -> x1(t) > 10, tt[t1_index+1:end])
    t2_index = t2_index !== nothing ? t2_index + t1_index : nothing
end

if t2_index !== nothing
    t3_index = findfirst(t -> x1(t) > 20, tt[t2_index+1:end])
    t3_index = t3_index !== nothing ? t3_index + t2_index : nothing
end

if t3_index !== nothing
    t4_index = findfirst(t -> x1(t) > 25, tt[t3_index+1:end])
    t4_index = t4_index !== nothing ? t4_index + t3_index : nothing
end

# Convert indices to times
t11 = t1_index !== nothing ? tt[t1_index] : "No such t1 found"
t22 = t2_index !== nothing ? tt[t2_index] : "No such t2 found"
t33 = t3_index !== nothing ? tt[t3_index] : "No such t3 found"
t44 = t4_index !== nothing ? tt[t4_index] : "No such t4 found"

println("First crossing time: ", t11)
println("Second crossing time: ", t22)
println("Third crossing time: ", t33)
println("Fourth crossing time: ", t44)
```

