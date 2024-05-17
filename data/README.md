# Data for replication

This folder contains the data needed for replicating the experiments of Section 6.1 (Conditional MNL instances). 

The instances of Section 6.2 (Generative MMNL instances) can be generated and added to this folder by running the script `scripts/generate_instances_MMNL.jl`. 

# File format:
The instances files have the following structure:

* Line 1: General parameters:  
<# customers (N)>, <# available facilities (D)>, <# competing facilities (E)>

* Line 2: Weight of each customer : 
<demand_customer_1> <demand_customer_2> ... <demand_customer_N>

* Line 3 to 2+E: Cost of each customer for the competing facilities: 
<cost_customer_1_comp_fac_1> ... <cost_customer_N_comp_fac_1> 
...
<cost_customer_1_comp_fac_E> ... <cost_customer_N_comp_fac_E> 

* Line 3+E to 2+E+D: Cost of each customer for the candidate locations: 
<cost_customer_1_cand_loc_1> ... <cost_customer_N_cand_loc_1> 
...
<cost_customer_1_cand_loc_D> ... <cost_customer_N_cand_loc_D> 
