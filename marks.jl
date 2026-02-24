
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





struct Schema <: Mark 
    μ :: Float64
    δ :: Float64
end

function ζ(sc::Schema)

    output = DiffusionRegion(μ, δ)
    output += Separatrix(μ, 1.0, :red)
    

end