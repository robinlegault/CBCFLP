# Author: Robin Legault <legault@mit.edu>

# Dependencies
include("../src/Experiments_runner.jl")

##########################################################################################################
################## Illustrative examples of Section 5 (Information-theoretic analysis) ###################
##########################################################################################################

# Parameters
inst_list = ["NYC_82341_59_5_1"]
β_list = [0.150, 0.145, 0.140, 0.135, 0.130, 0.125, 0.120, 0.115, 0.110, 0.105, 0.100, 0.095, 0.090, 0.085, 0.080, 0.075, 0.070, 0.065, 0.060, 0.055, 0.050]
α_list = [1]
S_list = [1]
b_list = [10]
excel_file = "entropy.xlsx"
std_file = "entropy.txt"
n_rep = 10
limit_time = 600
verbose = 2
dis = 2

# First call to CPLEX, avoid taking the overhead of the first call into account in the results
runner_MNL([inst_list[1]], [α_list[1]], [β_list[1]], [S_list[1]], [b_list[1]], "trash.xlsx"; verbose=0)

# Run experiments with Ω=Ω* (PBD method)
runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; tvd_min=0, clustering=true, P1_auto=true, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

# Run experiments with Ω=1 (Pure 0-1 method)
runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=false, P1_weight_min=1, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

# Run experiments with Ω=0 (Pure submodular method)
runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=false, P1_weight_min=0, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)