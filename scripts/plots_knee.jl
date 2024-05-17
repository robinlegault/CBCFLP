# Author: Robin Legault <legault@mit.edu>

# Packages
using DataFrames, XLSX, Statistics, Plots, LaTeXStrings, Random

# Dependencies
include("../src/IO_functions.jl")
include("../src/Simulation_and_clustering_functions.jl")

##########################################################################################################
############# Plots for illustrative examples of Section 4.2 (Partial Benders decomposition) #############
##########################################################################################################

# Load results
results_file = "knee_paper.xlsx"
df = DataFrame(XLSX.readtable(RESULTS_PATH*results_file, "Sheet1"; header=true, infer_eltypes=true))

# Results with predetermined relative weights Ω
df_pre = filter(row -> !row.P1_auto, df)
P1_weight_min_list = unique(df_pre.P1_weight_min)
n_pre = length(P1_weight_min_list)
mean_n_integer_nodes_BnC_pre = zeros(n_pre) # Number of integer solutions visited in the B&C    
mean_n_cuts_pre = zeros(n_pre)  # Number of submodular cuts generated
mean_time_pre = zeros(n_pre)    # Computing time (exact phase only, does not include the greedy heuristic)
mean_P1_pre = zeros(n_pre)      # Number of profiles in P1
for i ∈ 1:n_pre
    P1_weight_min = P1_weight_min_list[i]
    df_pre_i = filter(row -> row.P1_weight_min==P1_weight_min, df_pre)
    mean_n_integer_nodes_BnC_pre[i] = mean(df_pre_i.n_integer_nodes_BnC)
    mean_n_cuts_pre[i] = mean(df_pre_i.n_submodular_cuts)
    mean_time_pre[i] = mean(df_pre_i.time_exact)
    mean_P1_pre[i] = mean(df_pre_i.P1)
end

#Results with the knee method
df_auto = filter(row -> row.P1_auto, df)     
mean_n_integer_nodes_BnC_auto =  [mean(df_auto.n_integer_nodes_BnC)] # Number of integer solutions visited in the B&C    
mean_n_cuts_auto = [mean(df_auto.n_submodular_cuts)] # Number of submodular cuts generated
mean_time_auto = [mean(df_auto.time_exact)]          # Computing time (exact phase only, does not include the greedy heuristic)
mean_P1_auto = [mean(df_auto.P1)]                    # Number of profiles in P1
mean_P1_weight_auto = [mean(df_auto.P1_weight)]      # Relative weight Ω of P1

#Common attributes of the plots
default()
plot_font = "Computer Modern"
new_font = default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.2)
marker_size1 = 5
marker_size2 = 7
linewidth1 = 1
linewidth2 = 2
color_pre = "white"
color_auto = "black"
color_neutral1 = "grey"
color_neutral2 = "darkgrey"
dpi = 1000

#........................................................................................................#
#............ Value of δi by fraction of profiles included in P1 (illustrated on one instance) ..........#
#........................................................................................................#

# Load instance 1
Random.seed!(SEED)
α = 1
β = 0.1
S = 1
N, D, E, q, v, w_super, V, W_super = parse_MNL_instance("NYC_82341_59_5_1", α, β)
(a, ω) = simulation_MNL(v, w_super, q, S)
(a, ω) = exact_clustering(a, ω; remove_trivial=true)

# Compute P1 for each predetermined relative weights Ω of P1 
P = length(ω)
cum_weights = cumsum(ω)./sum(ω)
x_coord_list = zeros(length(P1_weight_min_list))
y_coord_list = zeros(length(P1_weight_min_list))
P1_list = zeros(Int, length(P1_weight_min_list))
for i in 1:length(P1_weight_min_list)
    P1_weight = P1_weight_min_list[i]
    if(P1_weight_min_list[i]==0)
        P1 = 0
    else
        P1 = findfirst(weight -> weight ≥ P1_weight, cum_weights)
        y_coord_list[i] = cum_weights[P1]
    end
    P1_list[i] = P1
    x_coord_list[i] = P1/P
end

# Limits of the plot
gap_plot_x = 0.03
x_lim_inf = 0 - gap_plot_x
x_lim_sup = 1 + gap_plot_x

# Plot points (i, Ωi) and (i, i)
x = (collect(0:P)./P)[1:100:end]
y = ([0] ∪ cum_weights)[1:100:end]
plt = plot(x, y, xlim = (x_lim_inf, x_lim_sup), linewidth=linewidth2, color=color_neutral1, label="Interpolated points "*L"(i, \Omega_i)", legend=:bottomright,dpi=dpi)
plot!(x, x, color=color_auto, label="Interpolated points "*L"(i, i)")
for i in 1:length(P1_weight_min_list)
    line_x = [x_coord_list[i], x_coord_list[i]]  # x-coordinates of the line segment
    line_y = [x_coord_list[i], y_coord_list[i]]  # y-coordinates of the line segment
    if(i==1)
        plot!(line_x, line_y, linestyle=:dash, color=color_neutral2, linewidth=linewidth2, label="Predetermined weights: "*L"δ_i")
    else
        plot!(line_x, line_y, linestyle=:dash, color=color_neutral2, linewidth=linewidth2)
    end
end

# Compute P1 using the knee method and plot (i*, Ωi*)
x = collect(range(1, P))./P
y = cum_weights
diffs = y-x
P1 = findlast(i -> i == maximum(diffs), diffs)
x_coord = P1/P
y_coord = cum_weights[P1]
line_x = [x_coord, x_coord]  # x-coordinates of the line segment
line_y = [x_coord, y_coord]  # y-coordinates of the line segment
plot!(line_x, line_y, linestyle=:dot, color=color_auto, linewidth=linewidth2, label="Knee: "*L"δ*")

# Axes labels
xlabel!(L"|\hat{P}_1| / \ |\hat{P}\ |")
ylabel!(L"Ω")

# x-axis ticks
custom_xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
custom_xticklabels = ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]
xticks!(custom_xticks, custom_xticklabels)

# Save the plot as a PNG file
savefig(FIGURES_OUTPUT_PATH*"plot_knee_choice.png")

#Print statistics
println("(|P|, i*, δ*, Ωi*) = (",P,", ",P1,", ",diffs[P1],", ",cum_weights[P1],")")


#........................................................................................................#
#............................... Computing time by relative weight Ω of P1 ..............................#
#........................................................................................................#

# x-axis values
x_pre = P1_weight_min_list 
x_auto = mean_P1_weight_auto

# y-axis values
y_pre = mean_time_pre
y_auto = mean_time_auto

# Insert point (x_auto, y_auto) in vectors x_pre and y_pre
pos_knee = findfirst(xi -> xi > x_auto[1], x_pre)
x_all = [x_pre[1:pos_knee-1]; x_auto[1]; x_pre[pos_knee:end]]
y_all = [y_pre[1:pos_knee-1]; y_auto[1]; y_pre[pos_knee:end]]

# Limits of the plot
gap_plot_y = 0.50
gap_plot_x = 0.03
y_lim_inf = 0
y_lim_sup = Int(ceil(maximum(y_all)))
x_lim_inf = 0 - gap_plot_x
x_lim_sup = 1 + gap_plot_x

# Horizontal line with average time of the knee method
plt = plot([x_lim_inf,x_lim_sup], [mean_time_auto[1],mean_time_auto[1]], linestyle=:dashdot, color=color_auto,
linewidth=linewidth2, label="", legend=:topleft, xlim = (x_lim_inf, x_lim_sup), ylim = (y_lim_inf, y_lim_sup), dpi=dpi)

# Line linking observed pairs (Ω, time)
plot!(x_all, y_all, label="", color=color_neutral1, linewidth=linewidth1)

# Observed pairs (Ω, time), predetermined weights
scatter!(x_pre, y_pre,  markersize=marker_size2, label="Predetermined weights", color=color_pre)

# Observed pair (average Ω, time), knee
scatter!(x_auto, y_auto, markersize=marker_size2, label="Knee", color=color_auto)

# Axes labels
xlabel!(L"Ω")
ylabel!("Average time "*L"(s)")

# x-axis and y-axis ticks
custom_xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
custom_xticklabels = ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]
xticks!(custom_xticks, custom_xticklabels)

custom_yticks = [round(y_auto[1], digits=2)]∪collect(y_lim_inf:y_lim_sup)
custom_yticklabels = string.(custom_yticks)
yticks!(custom_yticks, custom_yticklabels)

# Save the plot as a PNG file
savefig(FIGURES_OUTPUT_PATH*"knee_time.png")

#Print statistics
println("\nAverage time PBD = ",y_auto)
println("Average time pure 0-1 (Ω=1) = ",y_pre[end])
println("Average time pure submodular (Ω=0) = ",y_pre[1])


#........................................................................................................#
#........................ Number of submodular cuts and |P1| by relative weight Ω .......................#
#........................................................................................................#

# x-axis: Ω
x_pre = P1_weight_min_list 
x_auto = mean_P1_weight_auto

# y1-axis: average number of cuts, rescalesd for dual axes
scaling_factor = maximum(mean_P1_pre)/maximum(mean_n_cuts_pre)
y1_pre = mean_n_cuts_pre.*scaling_factor
y1_auto = mean_n_cuts_auto.*scaling_factor

# y2-axis: average cardinality of P1
y2_pre = mean_P1_pre
y2_auto = mean_P1_auto

# Insert points (x_auto, y_auto) in vectors x_pre and y_pre
pos_knee = findfirst(xi -> xi > x_auto[1], x_pre)
x_all = [x_pre[1:pos_knee-1]; x_auto[1]; x_pre[pos_knee:end]]
y1_all = [y1_pre[1:pos_knee-1]; y1_auto[1]; y1_pre[pos_knee:end]]
y2_all = [y2_pre[1:pos_knee-1]; y2_auto[1]; y2_pre[pos_knee:end]]

# Create the first plot with y1
plt = plot(x_all, y1_all, label="", color=color_neutral1, linewidth=linewidth1, right_margin=24Plots.mm, left_margin=2Plots.mm, size = (0.8*800, 0.8*500))
plot!(x_all, y2_all, label="", color=color_neutral1, linewidth=linewidth1, linestyle=:dash)
scatter!(x_pre, y1_pre, xlabel = L"Ω", ylabel = "Number of submodular cuts", label="Predetermined weights: #cuts", color=color_pre, legend=false, marker=:square, markersize=marker_size1,dpi=1000)
scatter!(x_auto, y1_auto, label="Knee: #cuts", color=color_auto, marker=:square, markersize=marker_size1)

# Add labels
scatter!(x_pre, y2_pre, ylabel = "Average number of cuts", label="Predetermined weights: "*L"| \hat{P}_1 \ |", markersize=marker_size2, color=color_pre)
scatter!(x_auto, y2_auto, label="Knee: "*L"| \hat{P}_1 \ |", secondary=true, markersize=marker_size2, color=color_auto)

# x-axis and y-axis ticks
custom_xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
custom_xticklabels = ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]
xticks!(custom_xticks, custom_xticklabels)
custom_yticks = [0,100,200,300,400].*scaling_factor
custom_yticklabels = ["0","100","200","300","400"]
yticks!(custom_yticks, custom_yticklabels)

# Create the second plot with y2 on a secondary y-axis
p = twinx()
scatter!(p, x_pre, y2_pre, ylabel = "Average cardinality of "*L"\hat{P}_1", label="", secondary=true, markersize=marker_size2, color=color_pre, yticks=[0,10000,20000,30000])
scatter!(p, x_auto, y2_auto, label="", secondary=true, markersize=marker_size2, color=color_auto)

# x-axis and y-axis ticks
custom_xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
custom_xticklabels = ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]
xticks!(custom_xticks, custom_xticklabels)

# Save the plot as a PNG file
savefig(FIGURES_OUTPUT_PATH*"knee_cuts.png")

#Print statistics
println("\nPBD: (#cuts, |P1|) = (",mean_n_cuts_auto,", ",mean_P1_auto,")")
println("PBD: # integer nodes B&C = ",mean_n_integer_nodes_BnC_auto)
println("Pure 0-1: (#cuts, |P1|) = (",mean_n_cuts_pre[end],", ",mean_P1_pre[end],")")
println("Pure submodular: (#cuts, |P1|) = (",mean_n_cuts_pre[1],", ",mean_P1_pre[1],")")
println("Pure submodular: # integer nodes B&C = ",mean_n_integer_nodes_BnC_pre[1])