
struct FieldLine <: Mark
    μ :: Float64
    Δ :: Float64
    r0 :: Float64
    w :: Float64
end
FieldLine(r0; μ=1.0, Δ=1.0, w=0.01) = FieldLine(μ, Δ, r0, w)

function ζ(fl::FieldLine)

    cy = (fl.r0 / fl.μ)^2
    pts = [(x, fl.μ*sqrt(x^2 + cy)) for x in LinRange(-fl.Δ/2, fl.Δ/2, 100)]
    return Trail(pts=pts, ws=fl.w)
end


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


struct Separatrix <: Mark
    μ :: Float64
    Δy :: Float64
    color :: Symbol
end
Separatrix(; μ=1.0, Δy=1.0, color=:black) = Separatrix(μ, Δy, color)

function ζ(sep::Separatrix)

    dx = sep.Δy/sep.μ
    ln1 = S(:strokeDasharray => 3, :fill => sep.color) * Line([-dx, dx], [-sep.Δy, sep.Δy])
    ln2 = S(:strokeDasharray => 3, :fill => sep.color) * Line([dx, -dx], [-sep.Δy, sep.Δy])

    return ln1 + ln2
end




struct DiffusionRegion <: Mark
    μ :: Float64
    δ :: Float64
end
DiffusionRegion(; μ=1.0, δ=1.0) = DiffusionRegion(μ, δ)

function ζ(dr::DiffusionRegion)

    rec = S(:opacity=>0.5, :fill=>:red)*Rectangle(h=dr.δ, w=dr.δ/dr.μ)
    arrow_1 = Arrow([0.0 ,0.0], [dr.δ, dr.δ/2], headsize=dr.δ/10)
    arrow_2 = Arrow([0.0 ,0.0], [-dr.δ, -dr.δ/2], headsize=dr.δ/10)
    arrow_3 = Arrow([dr.δ/dr.μ, (dr.δ+1)/dr.μ], [0.0, 0.0], headsize=dr.δ/10)
    arrow_4 = Arrow([-dr.δ/dr.μ, -(dr.δ+1)/dr.μ], [0.0, 0.0], headsize=dr.δ/10)

    return rec + arrow_1 + arrow_2 + arrow_3 + arrow_4
end


struct VerticalDist <: Mark
    Δx :: Float64
    Δy :: Float64
    label :: LaTeXString
    anchor :: Symbol
end
VerticalDist(Δy, label; anchor=:e, Δx=0.0) = VerticalDist(Δx, Δy, label, anchor)

function ζ(vd::VerticalDist)

    ln1 = S(:strokeDasharray => 7)*Line([0.0, 0.0], [vd.Δy/2, -vd.Δy/2])
    ln2 = S(:strokeDasharray => 3)*Line([-vd.Δx, 0.0], [vd.Δy/2, vd.Δy/2])
    ln3 = S(:strokeDasharray => 3)*Line([-vd.Δx, 0.0], [-vd.Δy/2, -vd.Δy/2])

    return ln1 + ln2 + ln3 + TextMark(text=vd.label, fontsize=vd.Δy/10, anchor=vd.anchor)
end


struct HorizontalDist <: Mark
    Δx :: Float64
    Δy :: Float64
    label :: LaTeXString
    anchor :: Symbol
end
HorizontalDist(Δx, label; anchor=:s, Δy=0.0) = HorizontalDist(Δx, Δy, label, anchor)

function ζ(hd::HorizontalDist)

    ln1 = S(:strokeDasharray => 7)*Line([-hd.Δx/2, hd.Δx/2], [0.0, 0.0])
    ln2 = S(:strokeDasharray => 3)*Line([-hd.Δx/2, -hd.Δx/2], [-hd.Δy, 0.0])
    ln3 = S(:strokeDasharray => 3)*Line([hd.Δx/2, hd.Δx/2], [-hd.Δy, 0.0])

    return ln1 + ln2 + ln3 + TextMark(text=hd.label, fontsize=hd.Δx/10, anchor=hd.anchor)
end


# ===============================
# ===============================

struct RGeodesic <: Mark
    K :: Float64
    d :: Float64
end

function ζ(rg::RGeodesic)

    sty = S(:strokeDasharray => 3)

    if rg.K == 0
        return sty*Line([-1.0, 1.0], [rg.d, rg.d])
    end

    r = 1/sqrt(abs(rg.K))
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

function ζ(gd::GeodesicGrid)


    if gd.K == 0.0
        crv = S(:strokeDasharray => 3)*Line([-1.0, 1.0], [0.0, 0.0])
    else
        crv = S(:strokeDasharray => 3)*T(0.0, -1.0)*Arc(abs(gd.K)^(-1/3), 1.0, [0.0, 0.0], 0.0, 0.0, π)
    end


    return crv
end








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
