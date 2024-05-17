# Author: Robin Legault <legault@mit.edu>

# Dependencies
include("../src/IO_functions.jl")
include("../src/Experiments_runner.jl")

##########################################################################################################
################## Computational experiments of Section 6.2 (Generative MMNL instances) ##################
##########################################################################################################

# Parameters
Dl = [25, 25, 25]
El = [10, 10, 10]
α  = 1
β_list = [[0.125, 0.25], [0.25, 0.5], [0.5, 1.0], [1.0, 2.0]]
N_list = [16000, 32000, 64000, 128000, 256000]
b_list = [5, 10, 15, 20, 25]
S_list = [5]
excel_file = "MMNL_MIX.xlsx"
std_file = "MMNL_MIX.txt"
limit_time = 600
verbose = 2
dis = 2
n_rep = 5

# Instance generation
for (β_min, β_max) ∈ β_list
    for N ∈ N_list
        for rep ∈ 1:n_rep
            instance = "MIX_N_"*string(N)*"_D1_"*string(Dl[1])*"_D2_"*string(Dl[2])*"_D3_"*string(Dl[3])*"_E1_"*string(El[1])*
                        "_E2_"*string(El[2])*"_E3_"*string(Dl[3])*"_α_"*string(α)*"_βmin_"*string(β_min)*"_βmax_"*string(β_max)*"_rep_"*string(rep)
            (N, D, E, q, v, w_super, V, W_super) = read_MNL_instance_utilities(instance)

            # Run experiments with PBD
            runner_MMNL(instance, q, v, V, w_super, W_super, b_list, S_list, excel_file; clustering=true, P1_auto=true, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

            # Run experiments with Ω=1 (Pure 0-1 method)
            runner_MMNL(instance, q, v, V, w_super, W_super, b_list, S_list, excel_file; clustering=true, P1_auto=false, P1_weight_min=1, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

            # Run experiments with SAA
            runner_MMNL(instance, q, v, V, w_super, W_super, b_list, S_list, excel_file; clustering=false, P1_auto=false, P1_weight_min=1, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)
        end
    end
end