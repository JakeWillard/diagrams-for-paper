# using Vizagrams
# import Vizagrams: ζ
using LaTeXStrings
using Plots

include("./formulas.jl")
include("./plots.jl")

plot_f_actual()



p = plot_lambda_0()
savefig(p, "./f_contours.png")
