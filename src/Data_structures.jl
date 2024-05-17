# Author: Robin Legault <legault@mit.edu>

##########################################################################################################
####################### High-level instance description and computational statistics #####################
##########################################################################################################
mutable struct Stats_PBD      # Description                                                # Default values 
    inst_name                 #Instance name                                               # ""
    D                         #Number of available locations                               # D (must be provided)
    b                         #Budget                                                      # 0
    N                         #Number of customers (conditional)                           # 0
    S                         #Number of scenarios (SAA)                                   # 0
    H                         #Entropy H(W_hat)                                            # 0.0
    Φ                         #Estimated error Φ_tilde(ω_hat)                              # 0.0
    δ                         #Optimal value δ_i* knee detection problem                   # 0.0
    weight_trivial            #Weight of the trivial preference pattern [0,0,...,0]        # 0.0
    weight_nontrivial         #Weight of the nontrivial preference patterns                # 0.0
    weight_total              #Weight of both trivial and nontrivial preference patterns   # 0.0
    time_clustering           #Clustering time (SAA -> DEQ_hat)                            # 0.0
    time_greedy               #Computing time greedy heuristic                             # 0.0
    time_exact                #Computing time exact phase                                  # 0.0
    n_nodes_BnC               #Visited nodes in the B&C tree                               # 0
    n_integer_nodes_BnC       #Visited integer nodes in the B&C tree                       # 0
    n_submodular_cuts         #Submodular cuts generated in the B&C tree                   # 0
    P                         #|P_hat|                                                     # 0
    P1                        #|P_hat_1|                                                   # 0
    P2                        #|P_hat_2|=|P_hat|-|P_hat_1|                                 # 0
    P1_auto                   #True is knee detection method used to select P1             # false
    P1_weight_min             #If P1_auto=false, minimum relative weight of P1             # 0.0
    P1_weight                 #(Weight P1)/(Weight P)                                      # 0.0
    x_greedy                  #Greedy solution                                             # []
    Z_greedy                  #Greedy value (LB)                                           # 0.0
    Z_greedy_stoch            #Objective value of x_greedy for stochastic model (e.g. MNL) # 0.0
    x_opt                     #Optimal solution                                            # []
    Z_opt                     #Optimal value                                               # 0.0
    Z_opt_stoch               #Objective value of x_opt for stochastic model (e.g. MNL)    # 0.0
    limit_time                #Time limit provided to the solver                           # 0.0
    is_optimal                #True if solved to proven optimality within the time limit   # false
end
function init_Stats_PBD(D ; inst_name="", b=0, N=0, S=0, H=0, Φ=0.0, δ=0.0, weight_trivial=0.0, weight_nontrivial=0.0,
     weight_total=0.0, time_clustering=0.0, limit_time=0.0)
     return Stats_PBD(inst_name, D, b, N, S, H, Φ, δ, weight_trivial, weight_nontrivial, weight_total,
     time_clustering, 0.0, 0.0, 0, 0, 0, 0, 0, 0,
     false, 0.0, 0.0, zeros(Bool, D), 0.0, 0.0, zeros(Bool, D), 0.0, 0.0, limit_time, false)
end

global column_names_PBD = [
"inst_name",                 
"D",                         
"b",                         
"N",                         
"S",                         
"H",                         
"Φ",                         
"δ",
"weight_trivial",
"weight_nontrivial",
"weight_total",      
"time_clustering",           
"time_greedy",               
"time_exact",                
"n_nodes_BnC",               
"n_integer_nodes_BnC",       
"n_submodular_cuts",         
"P",                         
"P1",                        
"P2",                        
"P1_auto",                   
"P1_weight_min",             
"P1_weight",                 
"x_greedy",                  
"Z_greedy",
"Z_greedy_stoch",                  
"x_opt",                     
"Z_opt",
"Z_opt_stoch",                     
"limit_time",                
"is_optimal"
]