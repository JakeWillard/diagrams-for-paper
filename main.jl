# using Vizagrams
# import Vizagrams: ζ
using LaTeXStrings
using Plots

include("./formulas.jl")
include("./plots.jl")


plot_orbit_vectors()


# p1 = light_rate_diagram(1.0, 1.0, (-2.0, 2.0))
# p2 = light_rate_diagram(0.5, 1.5, (-2.0, 2.0))
# p3 = light_rate_diagram(0.5, 0.5, (-2.0, 2.0))
# p4 = light_rate_diagram(1.5, 1.5, (-2.0, 2.0))

# title!(p1, "Special Relativistic")
# title!(p2, "Rightward Acceleration")
# title!(p3, "Tidal Compression")
# title!(p4, "Tidal Stretching")

# l = @layout [a b; c d]
# ptot = plot(p1, p2, p3, p4, layout=l, size=(800, 800))
# savefig(ptot, "./light-speed.png")


# light_rate_diagram(1.0, 1.0, (-2.0, 2.0))











# savefig(plot_rate(0.5), "./rates.png")

# savefig(plot_difference(), "./compare_rates.png")


# include("./marks.jl")





















# d = CurvScaleDiagram(1.5, 2.5, fs=0.2)
# d += T(10.0, 0.0) * U(0.5)*CurvScaleDiagram(2.5, 5.0, fs=0.4)
# draw(d)

# draw(CurvScaleDiagram(1.5, 2.5, fs=0.2))

# savefig(CurvScaleDiagram(1.5, 2.5, fs=0.2), filename="./bh-curvature.svg")



# tv = TensionVisual(0.0)
# tv += T(2.5, 0.0)*TensionVisual(0.4)
# tv += T(5, 0.0)*TensionVisual(-0.4)
# tv += T(0.0, 1.2)*TextMark(text=L"K=0", fontsize=0.2)
# tv += T(2.5, 1.2)*TextMark(text=L"K>0", fontsize=0.2)
# tv += T(5.0, 1.2)*TextMark(text=L"K<0", fontsize=0.2)

# draw(tv)

# savefig(tv, filename="./tension-illustration.svg")





# draw(FDDiagram(0.35))

# savefig(FDDiagram(0.35), filename="./schema.svg")

# savefig(RateIllustration(0.2), filename="./rate_illustration.svg")
# draw(RateIllustration(0.2))