# Author: Robin Legault <legault@mit.edu>

# Dependencies
include("../src/Instance_generation.jl")

##########################################################################################################
####### Generate instances for computational experiments of Section 6.2 (Generative MMNL instances) ######
##########################################################################################################

# Parameters
Dl = [25, 25, 25]
El = [10, 10, 10]
β_list = [[0.125, 0.25], [0.25, 0.5], [0.5, 1.0], [1.0, 2.0]]
α_list  = [1]
N_list = [125, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000, 1000000]
n_rep = 5

# Locations generation (same set of locations for all the instances)
Random.seed!(SEED)
coord_D, coord_E, type_D, type_E = generate_locations_mixture(Dl, El)

# Population generation (independent customers for all the instances)
for rep ∈ 1:n_rep
    for (β_min, β_max) ∈ β_list
        for α ∈ α_list
            for N ∈ N_list[1:end-1+(rep==1)] #only one rep for the instance of reference
                instance = "MIX_N_"*string(N)*"_D1_"*string(Dl[1])*"_D2_"*string(Dl[2])*"_D3_"*string(Dl[3])*"_E1_"*string(El[1])*
                                "_E2_"*string(El[2])*"_E3_"*string(Dl[3])*"_α_"*string(α)*"_βmin_"*string(β_min)*"_βmax_"*string(β_max)
                if(N < N_list[end])
                    instance = instance*"_rep_"*string(rep)
                end
                generate_instance_MMNL_mixture(N, coord_D, coord_E, type_D, type_E, α, instance, β_min, β_max)
            end
        end
    end
end