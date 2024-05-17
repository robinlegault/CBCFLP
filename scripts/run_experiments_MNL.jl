# Author: Robin Legault <legault@mit.edu>

# Dependencies
include("../src/Experiments_runner.jl")

##########################################################################################################
################## Computational experiments of Section 6.1 (Conditional MNL instances) ##################
##########################################################################################################

run_NYC = true
run_HM14 = true

if(run_NYC)
    # Parameters
    inst_list = ["NYC_82341_59_5_1"] 
    β_list = [0.10, 0.15, 0.20, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00]
    α_list = [0.75, 1, 1.25]
    S_list = [1, 5, 10]
    b_list = [2,3,4,5,6,7,8,9,10]
    excel_file = "MNL_NYC.xlsx"
    std_file = "MNL_NYC.txt"
    n_rep = 1
    limit_time = 3600
    verbose = 2
    dis = 2

    # First call to CPLEX, avoid taking the overhead of the first call into account in the results
    runner_MNL([inst_list[1]], [α_list[1]], [β_list[1]], [S_list[1]], [b_list[1]], "trash.xlsx"; verbose=0)

    # Run experiments with PBD
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=true, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

    # Run experiments with Ω=1 (Pure 0-1 method)
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=false, P1_weight_min=1, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

    # Run experiments with SAA
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=false, P1_auto=false, P1_weight_min=1, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)
end

if(run_HM14)
    # Parameters
    inst_list = ["HM14_800_25", "HM14_800_50", "HM14_800_100"] 
    β_list = [1, 2, 5, 10]
    α_list = [0.05, 0.1, 0.2]
    S_list = [10, 100, 1000]
    b_list = [2,3,4,5,6,7,8,9,10]
    excel_file = "MNL_HM14.xlsx"
    std_file = "MNL_HM14.txt"
    n_rep = 1
    limit_time = 600
    verbose = 2
    dis = 2

    # First call to CPLEX, avoid taking the overhead of the first call into account in the results
    runner_MNL([inst_list[1]], [α_list[1]], [β_list[1]], [S_list[1]], [b_list[1]], "trash.xlsx"; verbose=0)

    # Run experiments with PBD
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=true, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

    # Run experiments with Ω=1 (Pure 0-1 method)
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=true, P1_auto=false, P1_weight_min=1, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)

    # Run experiments with SAA
    runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; clustering=false, P1_auto=false, P1_weight_min=1, n_rep=n_rep, verbose=verbose, limit_time=limit_time, display=dis, std_file=std_file)
end