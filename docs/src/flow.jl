# Flow
using ForwardDiff, OrdinaryDiffEq

grad(f, x) = ForwardDiff.gradient(f, x)
jac(f, x) = ForwardDiff.jacobian(f, x)

function Flow(h)
    function rhs!(dz, z, dummy, t)
        n = size(z, 1) ÷ 2
        foo = z -> h(z[1:n], z[(n + 1):(2 * n)])
        dh = grad(foo, z)
        dz[1:n] = dh[(n + 1):(2n)]
        dz[(n + 1):(2n)] = -dh[1:n]
    end

    function f(tspan, x0, p0; abstol = 1e-12, reltol = 1e-12, saveat = [])
        z0 = [x0; p0]
        ode = ODEProblem(rhs!, z0, tspan)
        sol = OrdinaryDiffEq.solve(ode, Tsit5(), abstol = abstol, reltol = reltol, saveat = saveat)
        return sol
    end

    function f(t0, x0, p0, tf; abstol = 1e-12, reltol = 1e-12, saveat = [])
        sol = f((t0, tf), x0, p0, abstol = abstol, reltol = reltol, saveat = saveat)
        n = size(x0, 1)
        return sol[1:n, end], sol[(n + 1):(2 * n), end]
    end

    return f
end

function Lie(X, f)
    function Xf(x)
        df = grad(f, x)
        return df' * X(x)
    end

    return Xf
end

function Poisson(f, g)
    function fg(x, p)
        n = size(x, 1)
        ff = z -> f(z[1:n], z[(n + 1):(2n)])
        gg = z -> g(z[1:n], z[(n + 1):(2n)])
        df = grad(ff, [x; p])
        dg = grad(gg, [x; p])
        return df[(n + 1):(2n)]' * dg[1:n] - df[1:n]' * dg[(n + 1):(2n)]
    end

    return fg
end
