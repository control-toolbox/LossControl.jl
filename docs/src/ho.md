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
nothing # hide
```

```@example main
a0  = 0.0 
ε1  = 0.018
fNC(x) = fNC_unboundedminus(x,a0,ε1)
plot(fNC,-1., 1, label="fNC")
```

```@example main
@def ocp begin
        
    ε = 1e-3

    tf ∈ R,                          variable

    t ∈ [ 0., tf ],                  time

    q = [ x1, x2, λ ] ∈ R^3,         state

    ω = [u, v] ∈ R^2,                control

    tf ≥ 0.
    
    #initial conditions
    x1(0) == 2.5
    x2(0) == 4

    #final condition
    x1(tf) == 0
    x2(tf) == 0

    #control constraint
    -1. ≤  u(t)  ≤ 1.

    #state constraint
    -1 ≤  λ(t) ≤ 1,             (1)
    -5 ≤ x1(t) ≤ 5,             (2)
    -5 ≤ x2(t) ≤ 5,             (3)

    #control system
    q̇(t) == [x2(t), 
        (1-fNC(x2(t)))*u(t) + fNC(x2(t))*λ(t) - x1(t),
        (1-fNC(x2(t)))*v(t)]

    #cost function        
    tf + ∫(ε*(v(t))^2 +fNC(x2(t))*(u(t))^2) → min  
end
nothing # hide
```

```@example main
 N = 630 
sol = solve(ocp; init = (state = t -> [0.1, 0.1, 1, 0, 0], control =[-1, 0], variable =15), grid_size=N, print_level=4)
nothing # hide
```

```@example main
plot(sol; layout=:group, size=(800, 300))
```

```@example main
tf     = sol.variable
tt = (0:N+1) * value.(tf/(N+1))

x1(t)  = sol.state(t)[1]
x2(t)  = sol.state(t)[2]
λ(t)   = sol.state(t)[3]
u(t)   = sol.control(t)[1]
p1(t)  = sol.costate(t)[1]
p2(t)  = sol.costate(t)[2]
a      = λ(tf)
nothing # hide
```

```@example main
plot(x1, x2, 0, tf, label="optimal trajectory", color="blue", linewidth=2)
plot!([-4, 5], [0, 0], color=:black, label=false, linewidth=2)
```

```@example main
plot(tt, u, label="optimal control", color="red", linewidth=2)
plot!(tt, λ, label="state λ", color="green", linewidth=2)
```

```@example main
# Find the crossing times based on conditions for x1
t1_index = findfirst(t -> x2(t) ≤ 0, tt)
t2_index = nothing
t3_index = nothing

# If t1 is found, find the next crossing times
if t1_index !== nothing
    t2_index = findfirst(t -> x2(t) ≥ 0, tt[t1_index+1:end])
    t2_index = t2_index !== nothing ? t2_index + t1_index : nothing
end

if t2_index !== nothing
    t3_index = findfirst(t -> x2(t) ≤ 0, tt[t2_index+1:end])
    t3_index = t3_index !== nothing ? t3_index + t2_index : nothing
end

if t3_index !== nothing
    t4_index = findfirst(t -> x2(t) ≥ 0, tt[t3_index+1:end])
    t4_index = t4_index !== nothing ? t4_index + t3_index : nothing
end

# Convert indices to times
t1 = t1_index !== nothing ? tt[t1_index] : "No such t1 found"
t2 = t2_index !== nothing ? tt[t2_index] : "No such t2 found"
t3 = t3_index !== nothing ? tt[t3_index] : "No such t3 found"
t4 = tt[end]

println("First crossing time: ",        t1)
println("Second crossing time: ",       t2)
println("Third crossing time: ",        t3)
println("Fourth crossing/final time: ", t4)
```

```@example main
d = diff(u.(tt))
tstar = tt[findfirst(abs.(d) .> 0.9)[]]
println("the switching time: ", tstar)
```

```@example main
jmp1 = p2(t1+0.1)  - p2(t1-0.1)
jmp2 = p2(t2+0.1)  - p2(t2-0.1)

println(" p2(t1+) - p2(t1-) = ", jmp1)
println(" p2(t2+) - p2(t2-) = ", jmp2)
```

## Indirect Method 

```@example main
using NonlinearSolve  
using OrdinaryDiffEq
using Animations
nothing # hide
``` 

```@example main
# Dynamics
function F0(x)
    return [ x[2], -x[1]]
end

function F1(x)
    return [ 0.0 ,   1.0]
end

H0(x, p) = p' * F0(x) 
H1(x, p) = p' * F1(x)

# Hamiltonians: 
H(x, p, u)  =  H0(x, p) + u*H1(x,p)                             # pseudo-Hamiltonian

up(x, p) =   1.0
um(x, p) = - 1.0

Hp(x, p) = H(x, p, up(x, p))
Hm(x, p) = H(x, p, um(x, p))

# Hamiltonians: control loss region 2
H2(x, b, y, p)  = H0(x, p) + b*H1(x, p) - y*p[2]                # pseudo-Hamiltonian
Hcl(X, P)       = H2(X[1:2], X[3], X[4], P[1:2])                # control loss 2

# Flows
fp    = Flow(Hamiltonian(Hp))
fm    = Flow(Hamiltonian(Hm))
fcl   = Flow(Hamiltonian(Hcl))
nothing # hide
```

```@example main
# parameters
t0 = 0.0
x0 = [2.5; 4.0]
nothing # hide
```

```@example main 
# Shooting function
function shoot(p0, tt1, tt2, ttstar, tt3, b1, jump1, jump2, TT) 
    
    pb0    = 0.0 
    py0    = 0.0
        
    x1, p1 = fm(t0 , x0, p0, tt1) 
    
    x2, p2 = fp(tt1, x1, p1 - [0. , jump1], tt2)
    
    x3, p3 = fp(tt2, x2, p2, ttstar)
    
    x4, p4 = fm(ttstar, x3, p3, tt3)
    
    X5, P5 = fcl(tt3, [x4 ; b1 ; 0.0], [p4 - [0. , jump2]; pb0 ; py0], TT)

         s = zeros(eltype(p0), 10)
    
    s[1:2]  = X5[1:2] - [ 0.0 , 0.0 ]                     # target
    s[3]    = H1(x3, p3)                                  # switching
    s[4]    = x1[2] - 0.0                                 # first crossing 
    s[5]    = x2[2] - 0.0                                 # second crossing 
    s[6]    = x4[2] - 0.0                                 # third crossing 
    s[7]    = jump1 - p1[2]*(1. + 1.)/(1. - x1[1])        # jump1
    s[8]    = jump2 - p4[2]*(b1 + 1.)/(b1 - x4[1])        # jump2
    s[9]    = Hm(x0, p0) - 1.0                            # free final time
    s[10]   = P5[4]                                       # averaged gradient condition 

    return s

end
nothing # hide
``` 

```@example main
nle! =  (ξ, λ) -> shoot(ξ[1:2], ξ[3], ξ[4], ξ[5],ξ[6],ξ[7],ξ[8], ξ[9], ξ[10])
ξ_guess = [p1(0) , p2(0), t1, t2, tstar , t3, a, jmp1, jmp2, t4]; # initial guess
prob = NonlinearProblem(nle!, ξ_guess)
nothing # hide
```

```@example main
#solve
indirect_sol = solve(prob; abstol=1e-8, reltol=1e-8, show_trace=Val(true))
nothing # hide
```

```@example main
# Retrieves solution
pp0     = indirect_sol[1:2]
tt1     = indirect_sol[3]
tt2     = indirect_sol[4]
ttstar  = indirect_sol[5]
tt3     = indirect_sol[6]
b11     = indirect_sol[7]
jmp1    = indirect_sol[8]
jmp2    = indirect_sol[9]
T1      = indirect_sol[10]
nothing # hide
```

```@example main
ode_sol = fm((t0, tt1), x0, pp0, saveat=0.1) 
ttt1    = ode_sol.t 
xx1     = [ ode_sol[1:2, j] for j in 1:size(ttt1, 1) ] 
pp1     = [ ode_sol[3:4, j] for j in 1:size(ttt1, 1) ] 
uu1     = um.(xx1, pp1)

ode_sol = fp((tt1, tt2), xx1[end], pp1[end] - [0., jmp1], saveat=0.1) 
ttt2    = ode_sol.t ;
xx2     = [ ode_sol[1:2, j] for j in 1:size(ttt2, 1) ] 
pp2     = [ ode_sol[3:4, j] for j in 1:size(ttt2, 1) ] 
uu2     = up.(xx2, pp2)  

ode_sol = fp((tt2, ttstar), xx2[end], pp2[end] , saveat=0.1) 
ttt3    = ode_sol.t ;
xx3     = [ ode_sol[1:2, j] for j in 1:size(ttt3, 1) ] 
pp3     = [ ode_sol[3:4, j] for j in 1:size(ttt3, 1) ] 
uu3     = up.(xx3, pp3)  


ode_sol = fm((ttstar, tt3), xx3[end], pp3[end], saveat=0.1) 
ttt4    = ode_sol.t ;
xx4     = [ ode_sol[1:2, j] for j in 1:size(ttt4, 1) ] 
pp4     = [ ode_sol[3:4, j] for j in 1:size(ttt4, 1) ] 
uu4     = um.(xx4, pp4)  

ode_sol = fcl((tt3, T1), [xx4[end] ; b11 ; 0.0], [pp4[end] - [0., jmp2]; 0. ; 0.], saveat=0.1)
ttt5    = ode_sol.t
xx5     = [ ode_sol[1:2, j] for j in 1:size(ttt5, 1) ]
pp5     = [ ode_sol[5:6, j] for j in 1:size(ttt5, 1) ] 
uu5     = b11.*ones(length(ttt5)) 
nothing # hide
```

```@example main
tt0 = [ ttt1 ; ttt2 ; ttt3 ; ttt4 ; ttt5 ]
xx = [ xx1 ; xx2 ; xx3 ; xx4 ; xx5 ]
pp = [ pp1 ; pp2 ; pp3 ; pp4 ; pp5 ]
uu = [ uu1 ; uu2 ; uu3 ; uu4 ; uu5 ]

m = length(tt0)

x11 = [ xx[i][1] for i=1:m ]
x22 = [ xx[i][2] for i=1:m ]
p11 = [ pp[i][1] for i=1:m ]
p22 = [ pp[i][2] for i=1:m ]
nothing # hide
```

```@example main
plot(x11, x22, label="optimal trajectory",  linecolor=:blue , linewidth=2)
plot!([-4, 5], [0, 0], color=:black, label=false, linewidth=2)
```

```@example main
plot(tt0,   uu, label="optimal control",  linecolor=:red , linewidth=2) 
```

```@example main
plot(tt0,  p11, label="costate p1", linecolor=:purple , linewidth=2)
plot!(tt0,  p22, label="costate p2",  linecolor=:violet , linewidth=2)
```

```@example main
# create an animation
animx = @animate for i = 1:length(tt0)
    plot(x11[1:i], x22[1:i], xlim=(-3.,5.), ylim=(-4.,4.3), label="optimal trajectory", linecolor=:blue,  linewidth=2)
    scatter!([x11[i]], [x22[i]], markersize=4, marker=:circle, color=:black, label=false)
    plot!([-4, 5], [0, 0], color=:black, label=false, linewidth=2)

end
animu = @animate for i = 1:length(tt0)
    plot(tt0[1:i], uu[1:i], xlim=(0.,tt0[end]), ylim=(-1.2,1.2), label="opitmal control", linecolor=:red,  linewidth=2)
end 

animp1 = @animate for i = 1:length(tt0)
    plot(tt0[1:i], p11[1:i], xlim=(0.,tt0[end]), ylim=(-1.3, 0.5),label="costate p1", linecolor=:purple,  linewidth=2)
end ;

animp2 = @animate for i = 1:length(tt0)
    plot(tt0[1:i], p22[1:i], xlim=(0.,tt0[end]), ylim=(-1.5,1.3), label="costate p2", linecolor=:violet,  linewidth=2)
end
nothing # hide
```

```@example main
# display the animation
gif(animx, "ho_x.gif",   fps = 10)
```

```@example main
gif(animu, "ho_u.gif",   fps = 10)
```

```@example main
gif(animp2, "ho_p2.gif", fps = 10)
```
