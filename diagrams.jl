

function rate_illustration_diagram()


    V = 0.2
    x_stl = S(:strokeDasharray => 3, :stroke => :teal)
    y_stl = S(:strokeDasharray => 3, :stroke => :gold)

    geo = RateGeometry(V, x_stl=x_stl, y_stl=y_stl)

    Δt = abs((1 - sqrt(0.5^2 + 1))/V)
    h = Δt / 5
    brk1 = T(0.0, h) * Bracket(w=1.0, h=h)
    brk2 = T(0.75, -Δt/2) * R(-π/2) * Bracket(w=Δt, h=h)

    return draw(geo + brk1 + brk2)
end


function alpha_beta_definitions()

    ax = RotatedAxes(0.3)

    return draw(ax + T(1.5, 0.0) * ax)
end


function finite_difference_schema()

    fl_stl = S(:opacity => 0.2, :stroke => :darkblue)
    diff_stl = S(:opacity => 0.4, :fill => :grey)

    μ = 0.6
    δ = 0.1

    rec = Reconnection(μ, δ, fl_stl=fl_stl, diff_stl=diff_stl)
    
    sep = S(:strokeDasharray => 3) * Line([0.0, 1/μ], [0.0, 1.0])
    sep += S(:strokeDasharray => 3) * Line([0.0, -1/μ], [0.0, 1.0])

    shift = (1/μ)*1.5
    zstl = S(:strokeDasharray => 3, :stroke => :brown, :opacity => 0.5)
    zoom = ZoomIn(μ, δ, shift, zstl=zstl, dstl=diff_stl)

    r = 0.02
    fd_stl = S(:strokeDasharray => 7)
    upstream_cross = FDCross(μ, δ, 1.0, r, fd_stl)
    downstream_cross = T(shift, 0.0)*R(-π/2)*FDCross(1/μ, 0.0, 1/μ, r, fd_stl)

    sep += T(shift, 0.0)*S(:strokeDasharray => 3)*Line([0.0, 1/μ], [0.0, 1.0])
    sep += T(shift, 0.0)*S(:strokeDasharray => 3)*Line([0.0, 1/μ], [0.0, -1.0])

    bkts = T(0.0, 1.2)*Bracket(h=2*δ, w=1/μ)
    bkts += T(-0.75/μ, 0.5)*R(π/2)*Bracket(h=2*δ, w=1.0)

    return draw(rec + sep + zoom+ upstream_cross + downstream_cross + bkts)
end


function tension_illustration()

    grd_stl = S(:strokeDasharray => 2)
    ln_stl = S(:stroke => :darkblue)
    diff_stl = S(:opacity => 0.5, :fill => :grey)

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

    return draw(grds + lns + diff)
end


function light_diagram()



end


function f_visualizations()

    

end

