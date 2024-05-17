# Author: Robin Legault <legault@mit.edu>

# Dependencies
include("../src/Experiments_runner.jl")

##########################################################################################################
################## Illustrative examples of Section 4.2 (Partial Benders decomposition) ##################
##########################################################################################################

# Parameters
inst_list = ["NYC_82341_59_5_1"]
β_list = [0.1]
α_list = [1]
S_list = [1]
b_list = [10]
P1_weight_min_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
excel_file = "knee.xlsx"
std_file = "knee.txt"
n_rep = 10
verbose = 2
dis = 2

# First call to CPLEX, avoid taking the overhead of the first call into account in the results
runner_MNL([inst_list[1]], [α_list[1]], [β_list[1]], [S_list[1]], [b_list[1]], "trash.xlsx"; verbose=0)

# Run experiments with predetermined relative weights Ω of P1
for P1_weight_min ∈ P1_weight_min_list
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=false, P1_weight_min=P1_weight_min, n_rep=n_rep, verbose=verbose, display=dis, std_file=std_file)
end

# Run experiments with Ω=Ω* (knee method)
runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; tvd_min=0, clustering=true, P1_auto=true, n_rep=n_rep, verbose=verbose, display=dis, std_file=std_file)