using LaTeXStrings
using LinearAlgebra
# using ForwardDiff
rpath = "/home/jake/renders/"


# using Vizagrams
# import Vizagrams: ζ

using Plots

# using PlotlyJS

include("./formulas.jl")
# include("./marks.jl")
# include("./diagrams.jl")
include("./plots.jl")


savefig(plot_lambda(), rpath*"lambda.png")



plot_rate(0.3)
savefig(plot_rate(1.0), rpath*"rates.png")

savefig(plot_difference(), rpath*"thick_difference.png")




















only_tension_fieldlines("tension_only_fl_1.png", 0.005)
only_tension_fieldlines("tension_only_fl_2.png", 0.0075)
only_tension_fieldlines("tension_only_fl_3.png", 0.01)

only_fieldlines()




draw(Reconnection(0.6, 0.1))







savefig(plot_difference(), rpath * "diff_no_title.png")



savefig(plot_rate(0.5), rpath * "new_rates.png")

savefig(plot_lambda(), rpath * "lambdas.png")



new_tension_diagram()







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