# Author: Robin Legault <legault@mit.edu>

# Packages
using Random, Plots

# Dependencies
include("../src/Greedy_method.jl")
include("../src/IO_functions.jl")
include("../src/MNL_functions.jl")
include("../src/Plotting_functions.jl")
include("../src/PBD_method.jl")
include("../src/Simulation_and_clustering_functions.jl")

##########################################################################################################
########### Plots for illustrative examples of Section 3.2 (State-of-the-art heuristic method) ###########
##########################################################################################################

# Parameters
N  = 5000
S  = 1
α  = 0
α0 = 10
β  = 10
b  = 4

# Customers are distributed in the plane according to censored spherical normal centered at the origin
Random.seed!(SEED) 
coord_N = zeros(N,2)
μ =[0,0]
Σ = [[1.5 0.0];[0.0 1.5]]
for n=1:N
    n = Int(n)
    rep = true
    while (rep)
        coord_N[n,:] = rand(MvNormal(μ, Σ))
        if(coord_N[n,1]>=-2 && coord_N[n,1]<=2 && coord_N[n,2]>=-2 && coord_N[n,2]<=2)
            rep = false
        end
    end
end

# Location of the facilities
coord_E = [-1 1]
coord_D = [0 1 ; 1 1; -1 0; 0 0; 1 0; -1 -1; 0 -1; 1 -1]
D = size(coord_D)[1]
E = size(coord_E)[1]

# Generate instance
N, D, E, q, v, w_super, V, W_super = parse_MNL_instance_coord(coord_D, coord_E, coord_N, α, β; α0=α0)
(a, ω) = simulation_MNL(v, w_super, q, S)

# Greedy solution (suboptimal, identical to the solution returned by GGX; verified on the original MATLAB implementation)
(x_greedy, Z_greedy) = greedy_add(a, ω, b)
Z_greedy_MNL = objective_value_MNL(x_greedy,q,V,W_super)
println("x_greedy     = ",x_greedy)
println("Z_greedy     = ",Z_greedy,", (",round(100*Z_greedy/sum(ω), digits=2),"% of total)")
println("Z_greedy_MNL = ",Z_greedy_MNL,", (",round(100*Z_greedy_MNL/sum(ω), digits=2),"% of total)")

# PBD solution (optimal, identical to the solution returned by MOA; verified on the original MATLAB implementation)
(x_opt, Z_opt) = PBD(a, ω, b; return_stats=false)
Z_opt_MNL = objective_value_MNL(x_opt,q,V,W_super)
println("x_opt     = ",x_opt)
println("Z_opt     = ",Z_opt,", (",round(100*Z_opt/sum(ω), digits=2),"% of total)")
println("Z_opt_MNL = ",Z_opt_MNL,", (",round(100*Z_opt_MNL/sum(ω), digits=2),"% of total)")

# Generate plot for greedy solution and saves it as PNG files
plt_greedy = plot_solution_MNL(x_greedy, v, w_super, coord_N, coord_E, coord_D; legend=false)
savefig(FIGURES_OUTPUT_PATH*"GGX_example_suboptimal.png")

# Generate plot for optimal solution and saves it as PNG files
plt_opt = plot_solution_MNL(x_opt, v, w_super, coord_N, coord_E, coord_D; legend=false)
savefig(FIGURES_OUTPUT_PATH*"GGX_example_optimal.png")