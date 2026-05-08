using LaTeXStrings
using LinearAlgebra
using ForwardDiff
rpath = "/home/jake/renders/"


using Vizagrams
import Vizagrams: ζ

# using Plots

# using PlotlyJS

include("./formulas.jl")
include("./marks.jl")
include("./diagrams.jl")

only_geodesic_grid()
only_fieldlines()



# include("./rcm.jl")

conf = EMConfig((x,y) -> x^4 - y^4, 1.0)


Vf(0.0, 0.0, 0.01, conf)



cm(conf)

pts, l = exponential_map(0.0, 0.1, 0.01, conf)
l
plot_exponential_map(0.0, 0.4, 0.01, conf)




f_visualizations()





light_diagram()



tension_illustration()





finite_difference_schema()

alpha_beta_definitions()


rate_illustration_diagram()