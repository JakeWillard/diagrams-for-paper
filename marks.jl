
struct FieldLine <: Mark
    μ :: Float64
    Δ :: Float64
    r0 :: Float64
    w :: Float64
end
FieldLine(r0; μ=1.0, Δ=1.0, w=0.01) = FieldLine(μ, Δ, r0, w)

function ζ(fl::FieldLine)

    cy = (fl.r0 / fl.μ)^2
    Δx = sqrt((fl.Δ/fl.μ)^2 - cy)
    pts = [(x, fl.μ*sqrt(x^2 + cy)) for x in LinRange(-Δx, Δx, 100)]
    return Trail(pts=pts, ws=fl.w)
end


struct Separatrix <: Mark
    μ :: Float64
    Δy :: Float64
    color :: Symbol
end
Separatrix(; μ=1.0, Δy=1.0, color=:black) = Separatrix(μ, Δy, color)

function ζ(sep::Separatrix)

    dx = sep.Δy/sep.μ
    ln1 = S(:strokeDasharray => 3, :stroke => sep.color) * Line([-dx, dx], [-sep.Δy, sep.Δy])
    ln2 = S(:strokeDasharray => 3, :stroke => sep.color) * Line([dx, -dx], [-sep.Δy, sep.Δy])

    return ln1 + ln2
end


struct LabeledPoint <: Mark
    x :: Float64
    y :: Float64
    ϕ :: Float64
    label :: LaTeXString
    fontsize :: Float64
    dotsize :: Float64
end

function ζ(lp::LabeledPoint)

    pt = Circle(r=lp.dotsize, c=[lp.x, lp.y])
    dr = lp.dotsize + lp.fontsize/2.0
    txt = T(dr*cos(lp.ϕ), dr*sin(lp.ϕ))*T(lp.x, lp.y)*TextMark(text=lp.label, fontsize=lp.fontsize)
    return pt + txt
end


struct Upstream <: Mark
    μ :: Float64
    δ :: Float64
    Δy :: Float64
    fontsize :: Float64
end

function ζ(u::Upstream)

    dashes = S(:strokeDasharray => 7)
    greyout = S(:opacity => 0.5, :fill => :grey)

    Δy = u.Δy
    Δx = u.Δy / u.μ
    δ = u.δ
    l = δ / u.μ

    drgn = greyout*Rectangle(h=2*δ, w=2*l)

    pts = LabeledPoint(0.0, Δy/2, π/4, L"\mathfrak{q}", u.fontsize, u.fontsize/7)
    pts += LabeledPoint(0.0, Δy, π/4, L"\mathfrak{u}", u.fontsize, u.fontsize/7)
    pts += LabeledPoint(Δx/2, Δy/2, π/4, L"\mathfrak{h}", u.fontsize, u.fontsize/7)
    pts += LabeledPoint(0.0, δ, π/4, L"\mathfrak{m}", u.fontsize, u.fontsize/7)

    ylabel = dashes * Line([0.0, -3*Δx/4], [0.0, 0.0])
    ylabel += dashes * Line([0.0, 0.0], [0.0, Δy])
    ylabel += dashes * Line([0.0, -3*Δx/4], [Δy, Δy])
    ylabel += dashes * Line([-3*Δx/4, -3*Δx/4], [0.0, Δy])
    ylabel += T(-3*Δx/4, Δy/2) * TextMark(text=L"\Delta y^+", anchor=:w, fontsize=u.fontsize)

    xlabel = dashes * Line([-Δx/2, -Δx/2], [Δy/2, 3*Δy/2])
    xlabel += dashes * Line([-Δx/2, Δx/2], [Δy/2, Δy/2])
    xlabel += dashes * Line([Δx/2, Δx/2], [Δy/2, 3*Δy/2])
    xlabel += dashes * Line([-Δx/2, Δx/2], [3*Δy/2, 3*Δy/2])
    xlabel += T(0.0, 3*Δy/2) * TextMark(text=L"\Delta x^+", anchor=:s, fontsize=u.fontsize)


    sep = Separatrix(μ=u.μ, Δy=Δy, color=:black)

    nfl = 5
    r0s = LinRange(δ, Δy, nfl+1)[1:end-1]
    fls = greyout * FieldLine(r0s[1]; μ=u.μ, Δ=Δy, w=0.001)
    for i=2:nfl
        fls += greyout * FieldLine(r0s[i]; μ=u.μ, Δ=Δy, w=0.001)
    end
    fls += greyout * R(π)*FieldLine(r0s[1]; μ=u.μ, Δ=1.2*Δy, w=0.001)
    for i=2:nfl
        fls += greyout * R(π)*FieldLine((1 + 0.05*i)*r0s[i]; μ=u.μ, Δ=1.2*Δy, w=0.001)
    end
    #fls += R(π) * fls

    return drgn + pts + ylabel + xlabel + sep + fls
end

struct Downstream <: Mark
    μ :: Float64
    δ :: Float64
    fontsize :: Float64
end

function ζ(d::Downstream)

    Δy = 2*d.δ
    Δx = d.δ/d.μ

    drgn = S(:opacity => 0.3, :fill => :grey) * Rectangle(h=Δy, w=Δx)

    sep = S(:strokeDasharray => 3) * Line([-Δx/2, Δx/2], [0.0, Δy/2])
    sep += S(:strokeDasharray => 3) * Line([-Δx/2, Δx/2], [0.0, -Δy/2])

    pts = LabeledPoint(-Δx/2, 0.0, π, "x", d.fontsize, d.fontsize/7)
    pts += LabeledPoint(0.0, 0.0, 3*π/4, L"\mathfrak{p}", d.fontsize, d.fontsize/7)
    pts += LabeledPoint(0.0, Δy/4, 3*π/4, L"\mathfrak{s}", d.fontsize, d.fontsize/7)
    pts += LabeledPoint(Δx/2, 0.0, 0.0, L"\mathfrak{o}", d.fontsize, d.fontsize/7)

    dashes = S(:strokeDasharray => 7)

    xlabel = dashes * Line([-Δx/2, Δx/2], [0.0, 0.0])
    xlabel += dashes * Line([-Δx/2, -Δx/2], [0.0, 0.75*Δy])
    xlabel += dashes * Line([-Δx/2, Δx/2], [0.75*Δy, 0.75*Δy])
    xlabel += dashes * Line([Δx/2, Δx/2], [0.0, 0.75*Δy])
    xlabel += T(0.0, 0.75*Δy) * TextMark(text=L"l", anchor=:s, fontsize=d.fontsize)

    ylabel = dashes * Line([0.0, 0.0], [-Δy/4, Δy/4])
    ylabel += dashes * Line([0.0, 0.75*Δx], [Δy/4, Δy/4])
    ylabel += dashes * Line([0.0, 0.75*Δx], [-Δy/4, -Δy/4])
    ylabel += dashes * Line([0.75*Δx, 0.75*Δx], [-Δy/4, Δy/4])
    ylabel += T(0.75*Δx, 0.0) * TextMark(text=L"\delta", anchor=:e, fontsize=d.fontsize)


    return drgn + sep + pts + xlabel + ylabel
end


struct FDDiagram <: Mark
    μ :: Float64
    δ :: Float64
    Δy :: Float64
    fontsize :: Float64
end
FDDiagram(μ; δ=0.1, Δy=1.0, fontsize=0.15) = FDDiagram(μ, δ, Δy, fontsize)

function ζ(fd::FDDiagram)

    Δx = fd.Δy/fd.μ
    l = fd.δ/fd.μ
    xshift = 2Δx + l
    m =  (fd.Δy - fd.δ) / (xshift + Δx/2 - l)
    xw = xshift - Δx/2
    ynw = fd.Δy
    yne = fd.δ + m*(xshift - Δx/2 - l)
    ysw = -ynw
    yse = -yne

    u = Upstream(fd.μ, fd.δ, fd.Δy, fd.fontsize)
    d = T(xshift, 0.0)*Downstream(fd.μ, fd.Δy, fd.fontsize)

    lnstyle = S(:strokeDasharray => 3, :stroke => :purple)
    lns = lnstyle * Line([0.0, xw], [fd.δ, ynw])
    lns += lnstyle * Line([l, xw], [fd.δ, yne])
    lns += lnstyle * Line([0.0, xw], [-fd.δ, ysw])
    lns += lnstyle * Line([l, xw], [-fd.δ, yse])
    lns += lnstyle * Line([0.0, 0.0], [fd.δ, -fd.δ])
    lns += lnstyle * Line([l, l], [fd.δ, -fd.δ])
    lns += lnstyle * Line([0.0, l], [fd.δ, fd.δ])
    lns += lnstyle * Line([0.0, 0.0], [-fd.δ, -fd.δ])

    # xc = (fd.δ + 1.5*fd.Δy)/fd.μ
    # xw = (fd.δ + fd.Δy)/fd.μ
    # ynw = fd.Δy/2
    # ysw = -ynw

    return u + d + lns
end



# ===================================
# ===================================














# ===================================
# ===================================

struct FieldSheet <: Mark
    V :: Float64
    Δ :: Float64
    t0 :: Float64
    w :: Float64
end
FieldSheet(t0; V=1.0, Δ=1.0, w=0.01) = FieldSheet(V, Δ, t0, w)

function ζ(fs::FieldSheet)

    pts = [(s, (fs.V*fs.t0 + 1 - sqrt(s^2 + 1))/fs.V) for s in LinRange(-fs.Δ/2, fs.Δ/2, 100)]
    return Trail(pts=pts, ws=fs.w)
end


struct WorldLine <: Mark

    V :: Float64
    Δt :: Float64
    δ :: Float64
    t0 :: Float64
    w :: Float64
end
WorldLine(t0; V=1.0, Δt=1.0, δ=0.1, w=0.01) = WorldLine(V, Δt, δ, t0, w)

function ζ(wl::WorldLine)

    f(s) = 0.5*(tanh(s/wl.δ) + 1)
    t(s) = wl.t0 + (1 - f(s))*(wl.δ^2/wl.V)/s - f(s)*s/wl.V

    smax = wl.Δt * wl.V/2
    smin = wl.δ^2 / (wl.δ + smax)
    
    pts1 = [(s, t(s)) for s in LinRange(smin, smax, 100)]
    pts2 = [(-s, t(s)) for s in LinRange(smin, smax, 100)]

    return Trail(pts=pts1, ws=wl.w) + Trail(pts=pts2, ws=wl.w)
end


struct RateIllustration <: Mark

    V :: Float64

end

function ζ(r::RateIllustration)

    Δt = (sqrt(1.25)-1)/r.V

    fs = FieldSheet(0.0, V=r.V, w=0.003)
    fs += FieldSheet(Δt, V=r.V, w=0.003)
    fs += T(0.5, Δt/20)*TextMark(text=L"\psi_0", anchor=:se, fontsize=Δt/10)
    fs += T(0.5, Δt/20-Δt)*TextMark(text=L"\psi_0 \pm \Delta \psi", anchor=:se, fontsize=Δt/10)

    X = S(:stroke => :darkred, :strokeDasharray => 7) * Line([0.0, 0.0], [0.0, Δt])
    Y = S(:stroke => :darkgreen, :strokeDasharray => 7) * Line([-0.5, 0.5], [0.0, 0.0])

    clabels = S(:fill => :darkred)*T(0.0, Δt/2)*TextMark(text=L"X(\psi)", anchor=:w, fontsize=Δt/10)
    clabels += S(:fill => :darkgreen)*T(0.0, -Δt/20)*TextMark(text=L"Y", anchor=:n, fontsize=Δt/10)

    dlabels = S(:strokeDasharray => 3)*Line([0.0, 0.75], [Δt, Δt])
    dlabels += S(:strokeDasharray => 3)*Line([0.5, 0.75], [0.0, 0.0])
    dlabels += S(:strokeDasharray => 3)*Line([0.75, 0.75], [0.0, Δt])
    dlabels += T(0.75, Δt/2)*TextMark(text=L"w_t", anchor=:e, fontsize=Δt/10)
    dlabels += S(:strokeDasharray => 3)*Line([-0.5, -0.5], [0.0, -Δt/2])
    dlabels += S(:strokeDasharray => 3)*Line([0.5, 0.5], [0.0, -Δt/2])
    dlabels += S(:strokeDasharray => 3)*Line([-0.5, 0.5], [-Δt/2, -Δt/2])
    dlabels += T(0.0, -Δt/2)*TextMark(text=L"2w_y", anchor=:n, fontsize=Δt/10)

    return fs + X + Y + clabels + dlabels
end


# struct ConformalDiagram <: Mark
#     Cp :: Float64
#     Cm :: Float64
#     Δt :: Float64
#     Δ̄s :: Float64
#     fontsize :: Float64
# end
# ConformalDiagram(Cp, Cm; Δt=1.0, Δ̄s=1.0, fontsize=0.1) = ConformalDiagram(Cp, Cm, Δt, Δ̄s, fontsize)

# function ζ(cd::ConformalDiagram)

#     axes = Line([-cd.Δ̄s, cd.Δ̄s], [0.0, 0.0])
#     axes += Line([0.0, 0.0], [0.0, cd.Δt])
#     axes += T(cd.Δ̄s, 0.0)*TextMark(text="̂h", anchor=:ne, fontsize=cd.fontsize)
#     axes += T(0.0, cd.Δt)*TextMark(text=L"t", anchor=:nw, fontsize=cd.fontsize)

#     return axes
# end



# ===============================
# ===============================

struct RGeodesic <: Mark
    K :: Float64
    d :: Float64
end

function ζ(rg::RGeodesic)

    sty = S(:strokeDasharray => 3)

    if rg.K == 0 || rg.d == 0
        return sty*Line([-1.0, 1.0], [rg.d, rg.d])
    end

    r = 1/sqrt(rg.d*abs(rg.K))
    θ = atan(1/r)

    if rg.K > 0
        output = sty*Arc(r, r, [0.0, rg.d - r], 0.0, π/2-θ, π/2+θ)
    else
        output = sty*Arc(r, r, [0.0, r + rg.d], 0.0, 3*π/2-θ, 3*π/2+θ)
    end

    return output
end


struct GeodesicGrid <: Mark
    K :: Float64
    n :: Int64
end
GeodesicGrid(K; n=4) = GeodesicGrid(K, n)

function ζ(gg::GeodesicGrid)

    ds = LinRange(0.0, 1.0, gg.n)
    curves = RGeodesic(gg.K, ds[1])
    for i=2:gg.n
        curves += RGeodesic(gg.K, ds[i])
    end

    curves += R(π)*curves
    curves += R(π/2)*curves

    return curves
end


struct TensionVisual <: Mark
    K :: Float64
    δ :: Float64
    Δy :: Float64
    dcolor :: Symbol
    flcolor :: Symbol
end
TensionVisual(K; δ=0.1, Δy=0.05, dcolor=:grey, flcolor=:purple) = TensionVisual(K, δ, Δy, dcolor, flcolor)

function ζ(tv::TensionVisual)

    grd = GeodesicGrid(tv.K)
    drgn = S(:fill => tv.dcolor, :opacity => 0.5)*Rectangle(h=tv.δ, w=2.0)

    lns = S(:stroke => tv.flcolor) * Arrow([-1.0, 1.0], [tv.δ+tv.Δy, tv.δ+tv.Δy], headsize=tv.δ/7)
    lns += T(0.0, tv.Δy)*lns
    lns += T(0.0, 2*tv.Δy)*lns
    lns += R(π)*lns

    return grd + drgn + lns
end


# ================================
# ================================







# ===============================
# ===============================


struct SGeodesic <: Mark
    r :: Float64
    θ :: Float64
    w :: Float64
end

function ζ(sg::SGeodesic)

    Grrr(r) = 1/(2*r*(1-r))
    Grtt(r) = 1 - r
    Gtrt(r) = 1/r

    rs = sg.r*ones(2)
    ts = zeros(2)

    rdots = sqrt(1 - 1/sg.r)*cos(sg.θ)*ones(2)
    tdots = (sin(sg.θ)/sg.r)*ones(2)

    pts = []

    l = sg.r^(3/2) / 2.0
    ds = l/500
    for i=1:500
        rddots = -Grrr.(rs).*(rdots.^2) - Grtt.(rs).*(tdots.^2)
        tddots = -2*Gtrt.(rs).*rdots.*tdots
        rdots[2] = rdots[1] + ds*(rddots[1] + rddots[2])/2
        tdots[2] =  tdots[1] + ds*(tddots[1] + tddots[2])/2
        rs[2] = rs[1] + ds*(rdots[1] + rdots[2])/2
        ts[2] = ts[1] + ds*(tdots[1] + tdots[2])/2

        push!(pts, (rs[2]*cos(ts[2]), rs[2]*sin(ts[2])))

        if rs[2] <= 1.0
            break
        end

        rs[1] = rs[2]
        ts[1] = ts[2]
        rdots[1] = rdots[2]
        tdots[1] = tdots[2]
    end

    return Trail(pts=pts, ws=sg.w)
end


struct CurvatureVisual <: Mark
    r :: Float64
    ϕ :: Float64
    n :: Int64
    w :: Float64
    fs :: Float64
end
CurvatureVisual(r; ϕ=0.0, n=4, w=0.005, fs=0.1) = CurvatureVisual(r, ϕ, n, w, fs)

function ζ(cv::CurvatureVisual)


    ϵ = 1.5*cv.fs / (cv.r)
    output = S(:strokeDasharray => 7)*Arc(cv.r, cv.r, [0.0, 0.0], 0.0, -π/2+ϵ, 3*π/2-ϵ)

    output += T(0.0, -cv.r)*TextMark(text=L"%$(cv.r)r_s", fontsize=cv.fs)

    for θ in LinRange(0.0, 2*π, cv.n+1)[2:end]
        output += R(cv.ϕ) * SGeodesic(cv.r, θ, cv.w)
    end

    return output
end

struct CurvScaleDiagram <: Mark
    rmin :: Float64
    rmax :: Float64
    n :: Int64
    w :: Float64
    fs :: Float64
end
CurvScaleDiagram(rmin, rmax; n=4, w=0.01, fs=0.1) = CurvScaleDiagram(rmin, rmax, n, w, fs)

function ζ(csd::CurvScaleDiagram)

    rs = LinRange(csd.rmin, csd.rmax, 3)

    output = Circle()
    output += CurvatureVisual(rs[1], 0.0, csd.n, csd.w, csd.fs)
    output += CurvatureVisual(rs[2], π/2, csd.n, csd.w, csd.fs)
    output += CurvatureVisual(rs[3], π, csd.n, csd.w, csd.fs)

    return output
end
