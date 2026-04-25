using LaTeXStrings
using LinearAlgebra
rpath = "/home/jake/renders/"


using Vizagrams
import Vizagrams: ζ

# using Plots

# using PlotlyJS

include("./formulas.jl")
include("./marks.jl")
include("./diagrams.jl")

tension_illustration()





finite_difference_schema()

alpha_beta_definitions()


rate_illustration_diagram()