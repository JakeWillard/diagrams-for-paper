

# =====================================================

struct Curve <: Mark
    pts :: Vector{Tuple{Float64, Float64}}
    stl :: S
end
Curve(pts; stl=S()) = Curve(pts, stl)

function ζ(c::Curve)

    xs = [pt[1] for pt in c.pts]
    ys = [pt[2] for pt in c.pts]
    N = length(xs)

    output = c.stl * Line(xs[1:2], ys[1:2])
    for i=2:N-1
        output += c.stl * Line(xs[i:i+1], ys[i:i+1])
    end

    return output
end

struct Curve2 <: Mark
    pts :: Vector{Tuple{Float64, Float64}}
    ws :: Float64
    stl :: S
end
Curve2(pts; ws=0.1, stl=S()) = Curve2(pts, ws, stl)

function ζ(c::Curve2)

    return c.stl * Trail(c.pts, c.ws)
end

# =====================================================


struct Reconnection <: Mark
    μ :: Float64
    δ :: Float64
    fl_stl :: S
    diff_stl :: S
end
Reconnection(μ, δ; fl_stl=S(), diff_stl=S()) = Reconnection(μ, δ, fl_stl, diff_stl)

function fl_pts(m, y0, Δ)

    cy = (y0 / m)^2
    Δx = sqrt((Δ/m)^2 - cy)
    return [(x, m*sqrt(x^2 + cy)) for x in LinRange(-Δx, Δx, 100)]
end

function ζ(r::Reconnection)

    Δx = 1 / r.μ
    l = r.δ / r.μ

    diffusion_region = r.diff_stl * Rectangle(h=2*r.δ, w=2*l)

    y0s = LinRange(r.δ, 1.0, 5)[2:end]
    field_lines = Curve(fl_pts(r.μ, y0s[1], 1.0), stl=r.fl_stl)
    for y0 in y0s[2:end]
        field_lines += Curve(fl_pts(r.μ, y0, 1.0), stl=r.fl_stl)
    end
    for x0 in LinRange(l, 1/r.μ, 5)[2:end]
        field_lines += R(π/2)*Curve(fl_pts(1/r.μ, x0, 1.0/r.μ), stl=r.fl_stl)
    end
    field_lines += R(π) * field_lines


    return field_lines
end


# =====================================================


struct FDCross <: Mark

    m :: Float64
    ymin :: Float64
    Δ :: Float64
    r :: Float64
    stl :: S

end
FDCross(m, Δ, r; ymin=0.0, stl=S()) = FDCross(m, ymin, Δ, r, stl)

function ζ(c::FDCross)

    ym = c.Δ/2
    xm = ym / c.m

    output = Curve([(0.0, c.ymin), (0.0, c.Δ)], stl=c.stl)
    output += Curve([(-xm, ym), (xm, ym)], stl=c.stl)

    output += Circle(r=c.r, c=[0.0, c.ymin])
    output += Circle(r=c.r, c=[0.0, ym])
    output += Circle(r=c.r, c=[xm, ym])
    output += Circle(r=c.r, c=[-xm, ym])
    output += Circle(r=c.r, c=[0.0, c.Δ])

    return output
end

# =====================================================

struct ZoomIn <: Mark

    μ :: Float64
    δ :: Float64
    Δx :: Float64
    zstl :: S
    dstl :: S
end
ZoomIn(μ, δ, Δx; zstl=S(), dstl=S()) = ZoomIn(μ, δ, Δx, zstl, dstl)

function ζ(z::ZoomIn)

    l = z.δ/z.μ
    m1 = (1 - z.δ)/z.Δx
    m2 = (1 - z.δ)/(z.Δx + (1 - z.δ)/z.μ)
    y1 = z.δ + m1*z.Δx
    y2 = z.δ + m2*z.Δx


    output = Curve([(0.0, z.δ), (z.Δx, y1)], stl=z.zstl)
    output += Curve([(l, z.δ), (z.Δx, y2)], stl=z.zstl)
    output += Curve([(0.0, -z.δ), (z.Δx, -y1)], stl=z.zstl)
    output += Curve([(l, -z.δ), (z.Δx, -y2)], stl=z.zstl)
    output += Curve([(0.0, z.δ), (l, z.δ), (l, -z.δ), (0.0, -z.δ), (0.0, z.δ)], stl=z.zstl)
    output += T(z.Δx + 1.0/(2*z.μ), 0.0) * z.dstl * Rectangle(h=2.0, w=1.0/z.μ)

    return output
end


# =====================================================


struct FlowArrow <: Mark
    h :: Float64
    w :: Float64
    stl :: S
end
FlowArrow(; h=2.0, w=1.0, stl=S()) = FlowArrow(h, w, stl)

function ζ(a::FlowArrow)

    w = a.w
    h = a.h

    shape = Polygon([(w, 0.0), 
                    (0.5*w, 0.7*h), 
                    (1.0*w, 0.7*h), 
                    (0.0, h), 
                    (-w, 0.7*h),
                    (-0.5*w, 0.7*h),
                    (-w, 0.0)])

    return a.stl * shape
end


# =====================================================


struct InflowOutflow <: Mark

    μ :: Float64
    δ :: Float64
    fl_stl :: S
    diff_stl :: S
    in_stl :: S
    out_stl :: S

end
InflowOutflow(μ, δ; fl_stl=S(), diff_stl=S(), in_stl=S(), out_stl=S()) = InflowOutflow(μ, δ, fl_stl, diff_stl, in_stl, out_stl)

function ζ(f::InflowOutflow)

    rec = Reconnection(f.μ, f.δ, fl_stl=f.fl_stl, diff_stl=f.diff_stl)
    l = f.δ/f.μ

    arrows = T(0.0, 1+2*f.δ)*R(π)*FlowArrow(h=1.0, w=l, stl=f.in_stl)
    arrows += T(2*l, 0.0)*R(3*π/2)*FlowArrow(h=1.0, w=f.δ, stl=f.out_stl)
    arrows += R(π) * arrows
    

    return rec + arrows
end


# =====================================================


struct Bracket <: Mark
    h :: Float64
    w :: Float64
    stl :: S
end
Bracket(; h=1.0, w=1.0, stl=S()) = Bracket(h, w, stl)

function ζ(b::Bracket)

    xs = LinRange(-1, 1, 101)
    ys = [1 - tanh(10*abs(x)) - exp(10*(abs(x) - 1)) for x in xs]
    pts = [(b.w*xs[i]/2, b.h*ys[i]/2) for i=1:101]

    return Curve(pts, stl=b.stl)
end




# =====================================================



struct RateGeometry <: Mark
    V :: Float64
    fs_stl :: S 
    x_stl :: S
    y_stl :: S
end
RateGeometry(V; fs_stl=S(), x_stl=S(), y_stl=S()) = RateGeometry(V, fs_stl, x_stl, y_stl)


function ζ(g::RateGeometry)

    Δt = abs((1 - sqrt(0.5^2 + 1))/g.V)
    pts_1 = [(y, (1 - sqrt(y^2 + 1))/g.V) for y in LinRange(-0.6, 0.6, 100)]
    pts_2 = [(y, t-Δt) for (y,t) in pts_1[30:70]]

    field_sheets = Curve(pts_1, stl=g.fs_stl)
    field_sheets += Curve(pts_2, stl=g.fs_stl)

    X = Curve([(0.0, 0.0), (0.0, -Δt)], stl=g.x_stl)
    Y = Curve([(-0.5, -Δt), (0.5, -Δt)], stl=g.y_stl)

    return field_sheets + X + Y
end


# =====================================================


struct CoordinateAxes <: Mark
    hz :: Float64
    stl :: S
end
CoordinateAxes(; hz=0.05, stl=S()) = CoordinateAxes(hz, stl)

function ζ(a::CoordinateAxes)

    X = Arrow(pts=[(0.0, 0.0), (1.0, 0.0)], headsize=a.hz)
    Y = Arrow(pts=[(0.0, 0.0), (0.0, 1.0)], headsize=a.hz)

    return X + Y
end

# =====================================================

struct RotatedAxes <: Mark
    θ :: Float64
    r :: Float64
    hz :: Float64
    stl_1 :: S
    stl_2 :: S
    arc_stl :: S
end
RotatedAxes(θ; r=0.2, hz=0.05, stl_1=S(), stl_2=S(), arc_stl=S()) = RotatedAxes(θ, r, hz, stl_1, stl_2, arc_stl)

function ζ(a::RotatedAxes)

    output = CoordinateAxes(hz=a.hz, stl=a.stl_1)
    output += R(a.θ)*CoordinateAxes(hz=a.hz, stl=a.stl_2)

    output += Curve([(a.r*cos(t), a.r*sin(t)) for t in LinRange(0, a.θ, 100)], stl=a.arc_stl)
    output += Curve([(a.r*cos(t), a.r*sin(t)) for t in LinRange(π/2, π/2 + a.θ, 100)], stl=a.arc_stl)

    return output
end




# =====================================================


struct GeodesicGrid <: Mark
    K :: Float64
    N :: Int64
    stl :: S
end
GeodesicGrid(K; N=10, stl=S()) = GeodesicGrid(K, N, stl)

function geo_pts(K, y0)

    xmax = K < 0 ? sqrt(2(1 - abs(y0))/(abs(K)*abs(y0))) : 1.0
    xmax = xmax > 1.0 ? 1.0 : xmax
    return [(x,y0 - 0.5*K*y0*x^2) for x in LinRange(-xmax, xmax, 60)]
end

function ζ(g::GeodesicGrid)

    y0s = LinRange(-1, 1, g.N)
    grd = Curve(geo_pts(g.K, y0s[1]), stl=g.stl)
    for y0 in y0s[2:end]
        grd += Curve(geo_pts(g.K, y0), stl=g.stl)
    end
    grd += R(π/2) * grd

    return grd
end


# =====================================================

struct Dish <: Mark
    r1 :: Float64
    r2 :: Float64
    sld_stl :: S
    ln_stl :: S
end
Dish(r1, r2; sld_stl=S(), ln_stl=S()) = Dish(r1, r2, sld_stl, ln_stl)
Dish(r; sld_stl=S(), ln_stl=S()) = Dish(r/4, r, sld_stl, ln_stl)


function ζ(d::Dish)

    
    Δθ = π/3
    dsh_pts = [(d.r2*cos(t), d.r2*sin(t)) for t in LinRange(3*π/2-Δθ, 3*π/2+Δθ, 100)]
    push!(dsh_pts, dsh_pts[1])

    tip = d.sld_stl * Circle(r=d.r1)
    ant = d.ln_stl * Line([0.0, 0.0], [0.0, dsh_pts[1][2]])
    dsh = d.sld_stl * Polygon(dsh_pts)

    return tip + ant + dsh
end


# =====================================================

struct LightPaths <: Mark
    a :: Float64
    St :: Float64
    r :: Float64
    box_size :: Vector{Float64}
    ax_stl :: S
    dist_stl :: S
    dsh_stl :: S
    pth_stl :: S
end
LightPaths(a, St, r, box_size; ax_stl=S(), dist_stl=S(), dsh_stl=S(), pth_stl=S()) = LightPaths(a, St, r, box_size, ax_stl, dist_stl, dsh_stl, pth_stl)

function ζ(p::LightPaths)

    ybar(y) = y - 0.5*p.a*y^2 + (p.St - p.a^2)*y^3 /6
    ybp = ybar(1.0)
    ybm = ybar(-1.0)
    Δy, tmax = p.box_size

    margin = 1.2*p.r
    output = p.ax_stl * Line([-Δy/2-margin, -Δy/2-margin], [-margin, tmax])
    output += p.ax_stl * Line([-Δy/2-margin, Δy/2], [-margin, -margin])

    for yb in [ybar(y) for y in LinRange(-1.0, 1.0, 8)]
        output += p.dist_stl * Line([yb, yb], [-margin, tmax])
    end

    dsh = p.dsh_stl * Dish(p.r)
    output += dsh
    output += T(ybp, ybp) * R(π/4) * dsh
    output += T(ybm, -ybm) * R(-π/4) * dsh

    θ = abs(ybm) == abs(ybp) ? 0.0 : π/4
    output += T(0.0, 2*abs(ybm)) * R(π - θ) * dsh
    output += T(0.0, 2*abs(ybp)) * R(π + θ) * dsh

    output += Curve([(0.0, 0.0), (ybp, ybp), (0.0, 2*ybp)], stl=p.pth_stl)
    output += Curve([(0.0, 0.0), (ybm, -ybm), (0.0, -2*ybm)], stl=p.pth_stl)

    return output
end


# =====================================================


struct ConformalSpace <: Mark
    a :: Float64
    St :: Float64
    Δ̄y :: Float64
    N :: Int64
    axis_stl :: S
    dpts_stl :: S
end
ConformalSpace(a, St; Δ̄y=1.0, N=8, axis_stl=S(), dpts_stl=S()) = ConformalSpace(a, St, Δ̄y, N, axis_stl, dpts_stl)

function ζ(s::ConformalSpace)

    output = s.axis_stl * Line([0.0, 0.0], [0.0, 1.0])
    output += s.axis_stl * Line([-s.Δ̄y, s.Δ̄y], [0.0, 0.0])

    y = LinRange(-1.0, 1.0, s.N)
    ybars = y .- 0.5*s.a*y.^2 .+ (s.St - s.a^2)*y.^3 /6

    for ybar in ybars
        output += s.dpts_stl * Line([ybar, ybar], [0.0, 1.0])
    end

    return output
end

# =====================================================


struct LightSignals <: Mark
    a :: Float64
    St :: Float64
    stl :: S
end
LightSignals(a, St; stl=S()) = LightSignals(a, St, stl)

function ζ(s::LightSignals)

    ybar(y) = y - 0.5*s.a*y^2 + (s.St - s.a^2)*y^3 /6
    
    output = s.stl * Line([0.0, ybar(1.0)], [0.0, ybar(1.0)])
    output += s.stl * Line([ybar(1.0), 0.0], [ybar(1.0), 2*ybar(1.0)])
    output += s.stl * Line([0.0, ybar(-1.0)], [0.0, -ybar(-1.0)])
    output += s.stl * Line([ybar(-1.0), 0.0], [-ybar(-1.0), -2*ybar(-1.0)])

    return output
end


# =====================================================

struct Radii <: Mark
    rs :: Vector{Float64}
    Δx :: Float64
    stl :: S
end
Radii(rs, Δx; stl=S()) = Radii(rs, Δx, stl)
Radii(r::Float64, Δx; stl=S()) = Radii([r], Δx, stl)

function ζ(rad::Radii)

    output = Circle(r=1)
    for r in rad.rs
        θ = rad.Δx/r
        pts = [(r*cos(t), r*sin(t)) for t in LinRange(-π/2+θ/2, 3*π/2-θ/2, 100)]
        output += Curve(pts, stl=rad.stl)
    end

    return output
end