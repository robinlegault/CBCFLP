# Author: Robin Legault <legault@mit.edu>

#Packages
using Random

#Dependencies
include("Data_structures.jl")
include("Entropy_functions.jl")
include("IO_functions.jl")
include("MNL_functions.jl")
include("Simulation_and_clustering_functions.jl")
include("PBD_method.jl")

##########################################################################################################
################################ Runner for experiments on MNL instances #################################
##########################################################################################################
# inst_list    : list of instances, e.g. ["NYC_82341_59_5_1"]
# α_list       : list of aversion to competing facilities α
# β_list       : list of aversion to distances β
# S_list       : list of number of scenarios S
# b_list       : list of budgets b
# excel_file   : Excel file for collecting statistics, e.g. "results.xlsx"
#---------------------------------------------------------------------------------------------------------
# n_rep          : number of repetitions of the experiments, with different realizations of the r.v.s each time
# sheet_name     : append the row at the end of this sheet
# clustering     : if true, simulated customers with the same preference profile are aggregated (model \hat{DEQ}); otherwise, model SAA
# P1_auto        : if true, P1 is selected using the knee method (PBD)
# min_P_submod   : if |P|<min_P_submod, set P1=P. Applies if P1_auto is true
# tvd_min        : if δ*<tvd_min, set P1=P. Applies if P1_auto is true
# P1_weight_min  : if P1_auto=false, P1 is smallest set s.t. Ω ≥ P1_weight_min. 0->Pure submodular, 1->Pure 0-1
# heuristic      : if true, a greedy solution is first computed and provided to the solver
# verbose        : if 0, no standrd output on the instance. If 1, print instance id. If 2, print some results on the instance
# display        : if 0, no output is printed by the solver. Higher values provide more details
# limit_time     : time limit provided to the solver, in seconds
# std_file       : txt file on which the standard output is redirected, e.g. "results.xlsx"
##########################################################################################################
function runner_MNL(inst_list, α_list, β_list, S_list, b_list, excel_file; n_rep=1, sheet_name="Sheet1", clustering=true, 
    min_P_submod=0, tvd_min=0.0, P1_auto=true, P1_weight_min=1, heuristic=true, verbose=1, display=0, limit_time=3600, seed=SEED, std_file="")

    if(std_file!="")
        # Open a file in write mode
        stdfile = open(RESULTS_OUTPUT_PATH*std_file, "a")

        # Redirect stdout to the file
        redirect_stdout(stdfile)
    end

    #Number of instances to solve and counter
    n_inst = length(S_list)*length(inst_list)*length(α_list)*length(β_list)*length(b_list)*n_rep
    inst_count = 0

    # Instance
    for inst ∈ inst_list
        
        # Aversion to distances
        for β ∈ β_list

            # Aversion to competing facilities
            for α ∈ α_list

                # Load instance
                (N, D, E, q, v, w_super, V, W_super) = parse_MNL_instance(inst, α, β)

                # Number of scenarios
                for S ∈ S_list

                    #Common random numbers across all sets of parameters (S, α, β) to minimize variance
                    Random.seed!(seed) 

                    #Repetitions
                    for rep=1:n_rep
                    
                        # Draw a set of realization of the random terms
                        (a, ω) = simulation_MNL(v, w_super, q, S)
                        weight_total = sum(ω) # weight before removing the trivial preference profile
                        
                        # Exact clustering of the simulated customers. No clustering -> SAA
                        time_clustering = clustering ? (@timed begin (a, ω) = exact_clustering(a, ω; remove_trivial=true) end).time : 0
                        
                        # Budget
                        for b ∈ b_list
                            if(verbose≥1) println("__ Instance ",(inst_count+=1),"/",n_inst," __") end

                            # Initialize statistics
                            inst_name = inst*"_alpha_"*string(α)*"_beta_"*string(β)*"_S_"*string(S)*"_b_"*string(b)*"_rep_"*string(rep)
                            weight_nontrivial = sum(ω)
                            weight_trivial = weight_total - weight_nontrivial
                            H = empirical_entropy(ω)
                            Φ = empirical_estimated_expected_abs_dist(ω, sum(ω))
                            stats_PBD = init_Stats_PBD(D ; inst_name=inst_name, b=b, N=N, S=S, H=H, Φ=Φ, weight_trivial=weight_trivial, 
                            weight_nontrivial=weight_nontrivial, weight_total=weight_total, time_clustering=time_clustering, limit_time=limit_time)

                            # Solve the instance
                            stats_PBD = PBD(a, ω, b; P1_auto=P1_auto, min_P_submod=min_P_submod, tvd_min=tvd_min, P1_weight_min=P1_weight_min,
                             heuristic=heuristic, display=display, limit_time=limit_time, stats_PBD=stats_PBD)
                            
                            # Compute the objective value for the MNL problem
                            stats_PBD.Z_greedy_stoch = objective_value_MNL(stats_PBD.x_greedy,q,V,W_super)
                            stats_PBD.Z_opt_stoch = objective_value_MNL(stats_PBD.x_opt,q,V,W_super)

                            # Normalize the objective value of the simulation-based model according to the number of scenarios
                            stats_PBD.Z_greedy = stats_PBD.Z_greedy/S
                            stats_PBD.Z_opt    = stats_PBD.Z_opt/S

                            # Save the results in a file
                            save_results_PBD(stats_PBD, excel_file; sheet_name=sheet_name)

                            # Print basic information to standard output
                            if(verbose≥2) print_results_PBD(stats_PBD) end
                            flush(stdout)
                        end
                    end
                end
            end
        end
    end

    # Close the standard output file
    if(std_file!="")
        close(stdfile)
    end
end

##########################################################################################################
############################### Runner for experiments on MMNL instances #################################
##########################################################################################################
# instance     : instance name, e.g. "MIX_N_10000_D1_15_D2_15_D3_15_E1_5_E2_5_E3_5_α_1_βmin_1_βmax_2"
# q[n]         : weight of customer n
# V[n,d]       : exp(v[n,d]), where v[n,d] is the deterministic utility of location d for customer n
# W_super[n]   : exp(w_super[n]), where w_super[n] is the deterministic utility of the competing facility for customer n
# b_list       : list of budgets b
# S_list       : list of number of scenarios S
# excel_file   : Excel file for collecting statistics, e.g. "results.xlsx"
#---------------------------------------------------------------------------------------------------------
# sheet_name     : append the row at the end of this sheet
# clustering     : if true, simulated customers with the same preference profile are aggregated (model \hat{DEQ}); otherwise, model SAA
# P1_auto        : if true, P1 is selected using the knee method (PBD)
# min_P_submod   : if |P|<min_P_submod, set P1=P. Applies if P1_auto is true
# tvd_min        : if δ*<tvd_min, set P1=P. Applies if P1_auto is true
# P1_weight_min  : if P1_auto=false, P1 is smallest set s.t. Ω ≥ P1_weight_min. 0->Pure submodular, 1->Pure 0-1
# heuristic      : if true, a greedy solution is first computed and provided to the solver
# verbose        : if 0, no standrd output on the instance. If 1, print instance id. If 2, print some results on the instance
# display        : if 0, no output is printed by the solver. Higher values provide more details
# limit_time     : time limit provided to the solver, in seconds
# std_file       : txt file on which the standard output is redirected, e.g. "results.xlsx"
##########################################################################################################
function runner_MMNL(instance, q, v, V, w_super, W_super, b_list, S_list, excel_file; sheet_name="Sheet1", clustering=true, 
    min_P_submod=0, tvd_min=0.0, P1_auto=true, P1_weight_min=1, heuristic=true, verbose=1, display=0, limit_time=3600, seed=SEED, std_file="")

    N = length(q)
    D = length(v[1,:])

    if(std_file!="")
        # Open a file in write mode
        stdfile = open(RESULTS_OUTPUT_PATH*std_file, "a")

        # Redirect stdout to the file
        redirect_stdout(stdfile)
    end

    #Number of instances to solve and counter
    n_inst = length(b_list)
    inst_count = 0

    #Common random numbers across all sets of parameters (S, α, β) to minimize variance
    Random.seed!(seed) 

    # Number of scenarios
    for S ∈ S_list

        # Draw a set of realization of the random terms
        (a, ω) = simulation_MNL(v, w_super, q, S)
        weight_total = sum(ω) # weight before removing the trivial preference profile
        
        # Exact clustering of the simulated customers. No clustering -> SAA
        time_clustering = clustering ? (@timed begin (a, ω) = exact_clustering(a, ω; remove_trivial=true) end).time : 0
        
        # Budget
        for b ∈ b_list
            if(verbose≥1) println("__ Instance ",(inst_count+=1),"/",n_inst," __") end

            # Initialize statistics
            inst_name = instance*"_b_"*string(b)
            weight_nontrivial = sum(ω)
            weight_trivial = weight_total - weight_nontrivial
            H = empirical_entropy(ω)
            Φ = empirical_estimated_expected_abs_dist(ω, sum(ω))
            stats_PBD = init_Stats_PBD(D ; inst_name=inst_name, b=b, N=N, S=S, H=H, Φ=Φ, weight_trivial=weight_trivial, 
            weight_nontrivial=weight_nontrivial, weight_total=weight_total, time_clustering=time_clustering, limit_time=limit_time)

            # Solve the instance
            stats_PBD = PBD(a, ω, b; P1_auto=P1_auto, min_P_submod=min_P_submod, tvd_min=tvd_min, P1_weight_min=P1_weight_min,
                heuristic=heuristic, display=display, limit_time=limit_time, stats_PBD=stats_PBD)
            
            # Compute the objective value for the MNL problem
            stats_PBD.Z_greedy_stoch = objective_value_MNL(stats_PBD.x_greedy,q,V,W_super)
            stats_PBD.Z_opt_stoch = objective_value_MNL(stats_PBD.x_opt,q,V,W_super)

            # Normalize the objective value of the simulation-based model according to the number of scenarios
            stats_PBD.Z_greedy = stats_PBD.Z_greedy/S
            stats_PBD.Z_opt    = stats_PBD.Z_opt/S

            # Save the results in a file
            save_results_PBD(stats_PBD, excel_file; sheet_name=sheet_name)

            # Print basic information to standard output
            if(verbose≥2) print_results_PBD(stats_PBD) end
            flush(stdout)
        end
    end

    # Close the standard output file
    if(std_file!="")
        close(stdfile)
    end
end