
struct EMConfig

    Bx :: Function
    By :: Function
    Ez :: Function

end
function EMConfig(ψ::Function, E)

    Bx(x,y) = ForwardDiff.derivative(y -> ψ(x,y), y)
    By(x,y) = -ForwardDiff.derivative(x -> ψ(x,y), x)
    Ez(x,y) = E

    return EMConfig(Bx, By, Ez)
end

function g(x, y, conf::EMConfig)

    Bx = conf.Bx(x,y)
    By = conf.By(x,y)
    Ez = conf.Ez(x,y)

    gxx = 1 - By^2/Ez^2
    gxy = Bx*By/Ez^2
    gyy = 1 - Bx^2/Ez^2

    return [gxx gxy; gxy gyy]
end


function is_electric_dominant(x, y, conf::EMConfig)

    Cmin = 0.0
    C = conf.Ez(x,y)^2 - conf.Bx(x,y)^2 - conf.By(x,y)^2

    return C > Cmin
end


function Γ(x, y, conf::EMConfig)

    C1(x, y) = begin
        Bx = conf.Bx(x, y)
        By = conf.By(x, y)
        Ez = conf.Ez(x, y)
        return 0.5*[-By^2 0; 0 Bx^2] ./ Ez^2
    end

    C2(x,y) = begin
        Bx = conf.Bx(x, y)
        By = conf.By(x, y)
        Ez = conf.Ez(x, y)
        return 0.5*[0 -By^2; -By^2 2*Bx*By] ./ Ez^2
    end

    C3(x,y) = begin
        Bx = conf.Bx(x, y)
        By = conf.By(x, y)
        Ez = conf.Ez(x, y)
        return 0.5*[2*Bx*By -Bx^2; -Bx^2 0] ./Ez^2
    end

    Γx = ForwardDiff.derivative(x -> C1(x,y), x) + ForwardDiff.derivative(y -> C2(x,y), y)
    Γy = ForwardDiff.derivative(x -> C3(x,y), x) - ForwardDiff.derivative(y -> C1(x,y), y)

    return Γx, Γy
end

function f(vec, conf::EMConfig)

    x, y, xdot, ydot = vec

    if is_electric_dominant(x, y, conf::EMConfig) == false
        return [NaN, NaN, NaN, NaN]
    end

    Γx, Γy = Γ(x, y, conf)

    ginv = inv(g(x, y, conf))
    Γxup = ginv[1,1] * Γx + ginv[1,2] * Γy
    Γyup = ginv[2,1] * Γx + ginv[2,2] * Γy

    v = [xdot, ydot]
    xddot = -dot(v, Γxup * v)
    yddot = -dot(v, Γyup * v)

    return [xdot, ydot, xddot, yddot]
end


function rk4_step(vec, dt, conf::EMConfig)

    k1 = f(vec, conf)
    k2 = f(vec + dt*k1/2, conf)
    k3 = f(vec + dt*k2/2, conf)
    k4 = f(vec + dt*k3, conf)

    dv = dt*(k1 + 2*k2 + 2*k3 + k4)/6
    # if norm(dv[3:4]) / norm(vec[3:4]) > 1.0
    #     return [NaN, NaN, NaN, NaN]
    # end

    return vec + dv
end


function trace(vec0, dt, conf::EMConfig)

    pts = [(vec0[1], vec0[2])]
    vec = vec0[:]
    t = 0
    for dummy=1:1000
        vec[:] = rk4_step(vec[:], dt, conf)
        if isnan(vec[1])
            break
        end
        push!(pts, (vec[1], vec[2]))
        t += dt
    end

    return pts, t
end


function inits(x, y, conf::EMConfig)

    eigen_solve = eigen(g(x, y, conf))
    xhat = eigen_solve.vectors[:,1]
    yhat = eigen_solve.vectors[:,2]

    output = []
    for θ in LinRange(0, 2*π, 21)[1:20]

        v0 = cos(θ)*xhat + sin(θ)*yhat
        push!(output, [x, y, v0[1], v0[2]])
    end

    return output
end



function exponential_map(x, y, dt, conf::EMConfig)

    pts_vec = []
    dists = []

    for vec0 in inits(x, y, conf)

        pts, tf = trace(vec0, dt, conf)
        push!(pts_vec, pts)
        push!(dists, tf)
    end
    
    mean_dist = length(dists) > 0.0 ? sum(dists) / length(dists) : NaN
    return pts_vec, mean_dist
end


function plot_exponential_map(x, y, dt, conf::EMConfig)

    pts_vec, l = exponential_map(x, y, dt, conf)

    x = [pts[1] for pts in pts_vec[1]]
    y = [pts[2] for pts in pts_vec[1]]
    p = plot(x, y)
    for pts in pts_vec
        x = [v[1] for v in pts]
        y = [v[2] for v in pts]
        plot!(p, x, y)
    end

    return p
end

function Vf(x, y, dt, conf::EMConfig)

    N = 20

    eigen_solve = eigen(g(x, y, conf))
    xhat = eigen_solve.vectors[:,1]
    yhat = eigen_solve.vectors[:,2]

    v = [0.0, 0.0]

    for θ in LinRange(0, 2*π, N+1)[1:N]
        v0 = cos(θ)*xhat + sin(θ)*yhat
        pts, l = trace([x, y, v0[1], v0[2]], dt, conf)
        v += l*v0 / N
    end

    return v
end



function cm(conf::EMConfig)

    x = LinRange(-0.5, 0.5, 50)
    y = LinRange(-0.5, 0.5, 50)
    cms = zeros(50, 50)

    for i=1:50
        for j=1:50
            dummy, l = exponential_map(x[i], y[j], 0.01, conf)
            cms[j,i] = l
        end
        println(i)
    end

    return heatmap(x, y, cms)
end


