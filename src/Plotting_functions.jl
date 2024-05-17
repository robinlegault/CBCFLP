# Author: Robin Legault <legault@mit.edu>

# Packages
using DataFrames, XLSX, Statistics, Plots, LaTeXStrings, Random, Colors

# Dependencies
include("../src/IO_functions.jl")
include("../src/MNL_functions.jl")
include("../src/Simulation_and_clustering_functions.jl")

##########################################################################################################
############################## Visualization of a solution to a MNL instance #############################
##########################################################################################################
# x[d]    : 1 ⟺ location d is open in the solution
# v[n,d] : exp(v[n,d]), where v[n,d] is the deterministic utility of location d for customer n
# w[n]   : exp(w_super[n]), where w_super[n] is the deterministic utility of the competing facility for customer n
# coord_D : coordinates of candidate locations
# coord_E : coordinates of competing facilities
# coord_N : coordinates of customers
#---------------------------------------------------------------------------------------------------------
# legend           : if true, display the legend
# show_customers   : if true, display the location of the customers
# links_prob       : if true, display the probability for each customer to select each open facility
# color_prob       : if true, display the probability for each customer to select a location d∈D
# markersize_N     : size of the markers for customers
# markersize_D     : size of the markers for closed locations
# markersize_E     : size of the markers for competing facilities
# markersize_Dx    : size of the markers for open locations
# dpi              : resolution of the image
# color_D          : color of closed locations
# color_Dx         : color of open locations
# color_E          : color of competing facilities (for probability of selection)
##########################################################################################################
# plt : Plot
##########################################################################################################
function plot_solution_MNL(x, v, w, coord_N, coord_E, coord_D; 
    legend=true, show_customers=true, links_prob=false, color_prob=true, markersize_N=4, markersize_D=11, 
    markersize_E=17, markersize_Dx=15, dpi=1000, color_D="grey60", color_Dx="white", color_E="black")

    # Font
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=false)
    default()
    scalefontsizes(1.2)

    # Initialize the plot
    (N,D) = size(v)
    E = length(size(w)) == 1 ? 1 : size(w)[2]
    active = findall(xd->xd==1,x) 
    plt = scatter(dpi=dpi, size = (600, 500), xlab = L"\theta_1", ylabel = L"\theta_2", framestyle = :box)

    # Plot customers
    prob_MNL = probabilities_MNL(x,v,w)
    if(show_customers)
        # Illustrate the probabilites prob_MNL through links between customers and locations
        if(links_prob)
            for n in 1:N
                for e in 1:E
                    plot!([coord_N[n,1],coord_E[e,1]], [coord_N[n,2], coord_E[e,2]], color = color_E, label="", alpha=prob_MNL[n,D+e])
                end
                for d in 1:D
                    plot!([coord_N[n,1],coord_D[d,1]], [coord_N[n,2], coord_D[d,2]], color = color_Dx, label="", alpha=prob_MNL[n,d])
                end
            end
        end
        # Illustrate the probability of capturing each customer
        if(color_prob)
            scatter!(coord_N[:,1], coord_N[:,2], color=cgrad([color_E, color_Dx], [0, 0.5, 1]), zcolor=sum(prob_MNL[:,1:D],dims=2), markersize = markersize_N, label = "Customers")
        else
            scatter!(coord_N[:,1], coord_N[:,2], color = color_Dx, markersize = markersize_N, label = "Customers")
        end
    end

    # Plot the location of the facilities
    scatter!(coord_E[:,1], coord_E[:,2], color = color_E, markersize = markersize_E, shape =:diamond, label = "Competitor", markerstrokecolor="grey70", legend=legend, markerstrokewidth=3)
    scatter!(coord_D[:,1], coord_D[:,2], color = color_D, shape =:utriangle, markersize = markersize_D, label = "Closed locations", legend=legend, markerstrokewidth=3)
    scatter!(coord_D[active,1], coord_D[active,2], color = color_Dx, shape =:utriangle, markersize = markersize_Dx, label = "Open locations", legend=legend, markerstrokewidth=3)

    # Position of the legend
    if(legend) plot!(legend=:topleft) end

    # Reinitialize default font and return the plot
    default()
    return plt
end