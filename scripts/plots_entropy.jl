# Author: Robin Legault <legault@mit.edu>

# Packages
using DataFrames, XLSX, Statistics, Plots, LaTeXStrings, Random, StrRegex

# Dependencies
include("../src/IO_functions.jl")
include("../src/MNL_functions.jl")
include("../src/Simulation_and_clustering_functions.jl")

##########################################################################################################
############### Plots for illustrative examples Section 5 (Information-theoretic analysis) ###############
##########################################################################################################

# Load results
results_file = "entropy_paper.xlsx"
df = DataFrame(XLSX.readtable(RESULTS_PATH*results_file, "Sheet1"; header=true, infer_eltypes=true))
col_β = parse.(Float64, [match(r"beta_([\d.]+)", string).captures[1] for string in df.inst_name])
df[!, "β"] = col_β
list_β = reverse(sort!(unique(df.β)))
limit_time = df.limit_time[1]

# Collect statistics
df_knee   = filter(row -> row.P1_auto, df)
df_01     = filter(row -> (row.time_clustering>0.0 && !row.P1_auto && row.P1_weight_min==1) , df)
df_submod = filter(row -> (!row.P1_auto && row.P1_weight_min==0), df)
time_exact_knee   = zeros(length(list_β))
time_exact_01     = zeros(length(list_β))
time_exact_submod = zeros(length(list_β))
list_H_hat = zeros(length(list_β))
list_δ = zeros(length(list_β))
list_P = zeros(length(list_β))
list_weight_nontrivial = zeros(length(list_β))
for i ∈ 1:length(list_β)
    β = list_β[i]
    time_exact_knee[i]   = mean(filter(row -> row.β==β, df_knee).time_exact)
    time_exact_01[i]     = mean(filter(row -> row.β==β, df_01).time_exact)
    time_exact_submod[i] = mean(filter(row -> row.β==β, df_submod).time_exact)
    list_H_hat[i]        = mean(filter(row -> row.β==β, df_knee).H)
    list_δ[i]            = mean(filter(row -> row.β==β, df_knee).δ)
    list_P[i]            = mean(filter(row -> row.β==β, df_knee).P)
    list_weight_nontrivial[i] = mean(filter(row -> row.β==β, df_knee).weight_nontrivial)
end

# Instances solved within the time limit for each method
range_01 = findall(time -> time<limit_time, time_exact_01)
range_submod = findall(time -> time<limit_time, time_exact_submod)
range_knee = findall(time -> time<limit_time, time_exact_knee)

# Attributes of the plots
plot_font = "Computer Modern"
default()
new_font = default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.2)
marker_size1 = 5
marker_size2 = 7
linewidth1 = 1
linewidth2 = 2
color= "black"
linewidth = 3
dpi = 1000


#........................................................................................................#
#..................... Average computing time of each method for each level of entropy ..................#
#........................................................................................................#

# Ω=1 (pure 0-1 method)
plt = plot(list_H_hat[range_01], log.(time_exact_01[range_01]), label="Ω=1", linewidth=linewidth,
 linestyle=:solid, color=color, xlim=[7.9,10.424], ylim=log.([0.05,500]), legend=:topleft, dpi=dpi,right_margin=4Plots.mm)

# Ω=Ω* (knee method)
plot!(list_H_hat[range_knee], log.(time_exact_knee[range_knee]), label="Ω=Ω*", linewidth=linewidth,
linestyle=:dash, color=color)

# Ω=0 (pure submodular method)
plot!(list_H_hat[range_submod], log.(time_exact_submod[range_submod]), label="Ω=0",
 linewidth=linewidth, linestyle=:dot, color=color)

# Axes labels
xlabel!(L"H(\hat{W})")
ylabel!("Average time "*L"(s)")

# Customize the x-axis ticks with names
xtick_values = [8.0, 8.5, 9.0, 9.5, 10.0, 10.5]  # x-axis tick values
xtick_names = ["8.0", "8.5", "9.0", "9.5", "10.0", "10.5"]  # Corresponding tick names
xticks!(xtick_values, xtick_names)

# Customize the y-axis ticks with names
ytick_values = [log(0.1), log(1),log(10),log(100),log(500)]  # y-axis tick values
ytick_names = ["0.1","1","10","100","500"]  # Corresponding tick names
yticks!(ytick_values, ytick_names)

# Save the plot as a PNG file
savefig(FIGURES_OUTPUT_PATH*"plot_entropy.png")

#Print statistics
println("list_β                 = ",list_β)
println("list_H_hat             = ",list_H_hat)
println("list_δ                 = ",list_δ)
println("list_P                 = ",list_P)
println("list_weight_nontrivial = ",list_weight_nontrivial)
println("time_exact_01          = ",time_exact_01)
println("time_exact_submod      = ",time_exact_submod)
println("time_exact_knee        = ",time_exact_knee)


#........................................................................................................#
#............................... Solution quality for each level of entropy .............................#
#........................................................................................................#

# Load results MOA
results_file_MOA = "GGX_MOA_entropy_paper.xlsx"
df_MOA = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*results_file_MOA, "MOA"; header=true, infer_eltypes=true))

# Compute average MNL objective value of the optimal solution to \hat{DEQ}, S=1, for each value of β
instance = "NYC_82341_59_5_1"
α = 1
list_Z_SIM = zeros(length(list_β))
list_Z_MNL_SIM = zeros(length(list_β))
list_Z_MNL_MOA = zeros(length(list_β))
    
# For each level of randomness β
for i in 1:length(list_β)
    β = list_β[i]
    list_x_opt_β = []

    # Collect solution from MOA
    df_MOA_β = filter(row -> row.beta==β, df_MOA)
    str = df_MOA_β.solution[1]
    expr = Meta.parse(str)
    x_opt_MOA = eval(expr)
    
    # Collect solutions from simulation-based method, S=1
    for str in filter(row -> row.β==β, df_01).x_opt 
        expr = Meta.parse(str)
        x_opt = eval(expr)
        push!(list_x_opt_β, x_opt)
    end
    
    # Load MNL instance
    N, D, E, q, v, w_super, V, W_super = parse_MNL_instance(instance, α, β)
    
    # Report the average optimal value of model \hat{DEQ} on this set of instances
    list_Z_SIM[i] = mean(filter(row -> row.β==β, df_01).Z_opt)

    # Compute the average value, for model DEQ, of the optimal solutions to model \hat{DEQ} on this set of instances
    list_Z_MNL_SIM[i] = mean([objective_value_MNL(x_opt,q,V,W_super) for x_opt ∈ list_x_opt_β])

    # Compute the optimal value of model DEQ
    list_Z_MNL_MOA[i] = objective_value_MNL(x_opt_MOA,q,V,W_super)

    println("")
    println("β       = ",β)
    println("V[1,:] = ",V[1,:])
    println("x_opt_MOA = ",x_opt_MOA)
    println("list_Z_MNL_SIM[i] = ",list_Z_MNL_SIM[i])
    println("list_Z_MNL_MOA[i] = ",list_Z_MNL_MOA[i])
end

#Print statistics
println("list_Z_MNL_SIM = ",list_Z_MNL_SIM)
println("list_Z_MNL_MOA = ",list_Z_MNL_MOA)
println("Z_SIM      = ",list_Z_SIM)
println("Z_MNL_SIM  = ",list_Z_MNL_SIM)
println("Z_MNL_MOA = ",list_Z_MNL_MOA)
println("Optimality gap (%) = ",((list_Z_MNL_MOA-list_Z_MNL_SIM)./list_Z_MNL_MOA).*100)