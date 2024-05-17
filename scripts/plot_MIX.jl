# Author: Robin Legault <legault@mit.edu>

# Packages
using Random, Plots

# Dependencies
include("../src/Plotting_functions.jl")
include("../src/Instance_generation.jl")

##########################################################################################################
################# Plot for MMNL-MIX instances of Section 6.2 (Generative MMNL isntances) #################
##########################################################################################################

# Parameters
Dl = [25, 25, 25]
El = [10, 10, 10]
N = 10000

# Locations generation (same set of locations for all the instances)
Random.seed!(SEED)
coord_D, coord_E, type_D, type_E = generate_locations_mixture(Dl, El)

D = length(type_D)
E = length(type_E)
coord_N, type_N, q = generate_population_mixture(N)

color_D  = "white"
color_E  = "black"
colors_N1 = ["grey20", "grey50", "grey90"]
markershapes = [:rect, :diamond, :hexagon, :star5, :xcross, :circle, :cross, :utriangle, :dtriangle, :rtriangle, :ltriangle,  :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

# Font
plot_font = "Computer Modern"
default()
default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.2)

# Scatter plot
plt = scatter(coord_N[:,1], coord_N[:,2], color=colors_N1[type_N], markersize = 2, markershape=:circle, label = "", dpi=1000, size = (600, 500), xlab = L"\theta_1", ylabel = L"\theta_2", ylims=[-17.5,16], xlims=[-17.5,16], framestyle = :box, legend = :outerbottom, legendcolumns = 4)
scatter!(coord_E[:,1], coord_E[:,2], color = color_E, markersize = 8, markershape=markershapes[type_E], label = "", markerstrokecolor="grey70", markerstrokewidth=2)
scatter!(coord_D[:,1], coord_D[:,2], color = color_D, markersize = 8, markershape=markershapes[type_D], label = "", markerstrokewidth=2)

# Saves the plot
savefig(FIGURES_OUTPUT_PATH*"MMNL_MIX_"*string(N)*".png")