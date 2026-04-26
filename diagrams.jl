

function rate_illustration_diagram()


    V = 0.2
    x_stl = S(:strokeDasharray => 3, :stroke => :darkred)
    y_stl = S(:strokeDasharray => 3, :stroke => :darkgreen)

    geo = RateGeometry(V, x_stl=x_stl, y_stl=y_stl)

    Δt = abs((1 - sqrt(0.5^2 + 1))/V)
    h = Δt / 5
    brk1 = T(0.0, h) * Bracket(w=1.0, h=h)
    brk2 = T(0.75, -Δt/2) * R(-π/2) * Bracket(w=Δt, h=h)

    output = geo + brk1 + brk2
    savefig(output, filename="rate_illustration.svg", directory=rpath)
    return draw(geo + brk1 + brk2)
end


function alpha_beta_definitions()

    ax = RotatedAxes(0.3)

    output = ax + T(1.5, 0.0) * ax
    savefig(output, filename="alpha_beta.svg", directory=rpath)
    return draw(output)
end


function finite_difference_schema()

    fl_stl = S(:opacity => 0.2, :stroke => :darkblue)
    diff_stl = S(:opacity => 0.4, :fill => :orange)

    μ = 0.6
    δ = 0.1

    rec = Reconnection(μ, δ, fl_stl=fl_stl, diff_stl=diff_stl)
    
    sep = S(:strokeDasharray => 3) * Line([0.0, 1/μ], [0.0, 1.0])
    sep += S(:strokeDasharray => 3) * Line([0.0, -1/μ], [0.0, 1.0])

    shift = (1/μ)*1.5
    zstl = S(:strokeDasharray => 3, :stroke => :orange, :opacity => 0.5)
    zoom = ZoomIn(μ, δ, shift, zstl=zstl, dstl=diff_stl)

    r = 0.03
    fd_stl = S(:strokeDasharray => 2)
    upstream_cross = FDCross(μ, δ, 1.0, r, fd_stl)
    downstream_cross = T(shift, 0.0)*R(-π/2)*FDCross(1/μ, 0.0, 1/μ, r, fd_stl)

    sep += T(shift, 0.0)*S(:strokeDasharray => 3)*Line([0.0, 1/μ], [0.0, 1.0])
    sep += T(shift, 0.0)*S(:strokeDasharray => 3)*Line([0.0, 1/μ], [0.0, -1.0])

    bkts = T(0.0, 1.2)*Bracket(h=2*δ, w=1/μ)
    bkts += T(-0.75/μ, 0.5)*R(π/2)*Bracket(h=2*δ, w=1.0)

    output = rec + sep + zoom+ upstream_cross + downstream_cross + bkts
    savefig(output, filename="shema.svg", directory=rpath)
    return draw(output)
end


function tension_illustration()

    grd_stl = S(:strokeDasharray => 1, :opacity=>0.3)
    ln_stl = S(:stroke => :blue)
    diff_stl = S(:opacity => 0.3, :fill => :orange)

    grds = GeodesicGrid(0.0, stl=grd_stl)
    grds += T(2.4, 0.0)*GeodesicGrid(0.7, stl=grd_stl)
    grds += T(4.8, 0.0)*GeodesicGrid(-2.0, stl=grd_stl)

    δ = 0.1
    lns = ln_stl * Line([-1.0, 1.0], [2*δ, 2*δ])
    for n=3:6
        lns += ln_stl * Line([-1.0, 1.0], [n*δ, n*δ])
    end
    lns += R(π)*lns
    lns += T(2.4, 0.0)*lns + T(4.8, 0.0)*lns

    diff = diff_stl * Rectangle(h=δ, w=2.0)
    diff += T(2.4, 0.0)*diff + T(4.8, 0.0)*diff

    output = grds + lns + diff
    savefig(output, filename="tension_illustration.svg", directory=rpath)
    return draw(output)
end


function light_diagram()

    Δx = 5.5
    Δy = 5.5
    box_size = [3.5, 3.5]

    r = 0.25
    dist_stl = S(:strokeDasharray => 3, :opacity => 0.5)
    pth_stl = S(:strokeDasharray => 7, :stroke => :darkred)

    dgm1 = LightPaths(0.0, 0.0, r, box_size, dist_stl=dist_stl, pth_stl=pth_stl)
    dgm2 = LightPaths(0.5, 0.0, r, box_size, dist_stl=dist_stl, pth_stl=pth_stl)
    dgm3 = LightPaths(0.0, 4.0, r, box_size, dist_stl=dist_stl, pth_stl=pth_stl)
    dgm4 = LightPaths(0.0, -2.0, r, box_size, dist_stl=dist_stl, pth_stl=pth_stl)

    output = T(0.0, Δy)*dgm1
    output += T(Δx, Δy)*dgm2
    output += T(0.0, 0.0)*dgm3
    output += T(Δx, 0.0)*dgm4

    savefig(output, filename="light_paths.svg", directory=rpath)
    return draw(output)
end


function f_visualizations()

    r1 = 5.0
    r2 = 7.0
    fl_stl = S(:stroke => :darkblue)
    diff_stl = S(:fill => :orange)
    in_stl = S(:fill => :darkgreen)
    out_stl = S(:fill => :darkred)
    hz = 0.1

    bh = Circle(r=1.0)

    lns = S(:strokeDasharray => 3) * Line([0.0, r2], [0.0, 0.0])
    lns += R(π/3) * lns
    lns += R(2*π/3) * S(:strokeDasharray => 3) * Line([0.0, r1], [0.0, 0.0])

    rec_shp = InflowOutflow(0.6, 0.1, fl_stl=fl_stl, diff_stl=diff_stl, in_stl=in_stl, out_stl=out_stl)
    recs = T(r1, 0.0)*rec_shp
    recs += R(π/3) * T(r1, 0.0) * R(π/2) * rec_shp

    ax = T(r2, 0.0)*CoordinateAxes(hz=hz)
    ax += R(π/3) * ax

    plane = T(r1, 0.0) * S(:fill => :darkblue) * Rectangle(h=2.0, w=0.2)
    plane += Arrow(pts=[(r1, 0.0), (r1+1.0, 0.0)], headsize=hz)
    plane = R(2*π/3) * plane

    output = bh + lns + recs + ax + plane
    savefig(output, filename="f_illustration.svg", directory=rpath)
    return draw(output)
end


function scale_visualization()

    κ = 0.1
    rs = [2.0, 4.0, 6.0]
    θs = LinRange(0, π, 3)
    scales = [sqrt(κ*r^3) for r in rs]

    rstl = S()
    fl_stl = S()
    diff_stl = S()

    output = Radii(rs, 2.0, stl=rstl)
    rec = R(π/2)*Reconnection(0.6, 0.1, fl_stl=fl_stl, diff_stl=diff_stl)
    for i=1:3
        r = rs[i]
        θ = θs[i]
        s = scales[i]
        output += R(θ) * T(r, 0.0) * U(s) * rec
    end

    return draw(output)
end

