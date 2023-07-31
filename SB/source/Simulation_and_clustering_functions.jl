# Author: Robin Legault <legault@mit.edu>

#Packages
using Distributions, Random

#Dependencies
include("Constants.jl")

##########################################################################################################
################################## Simulation for conditional MLN model ##################################
############################### Generate model SAA(NS) if clustering=false ###############################
############################### Generate model \hat{DEQ} if clustering=true ##############################
##########################################################################################################
# v[n,d] : deterministic utility of location d for customer n
# w[n,c] : deterministic utility of competitor c for customer n
# q[n]   : weight of customer n
# S      : number of scenarios
#---------------------------------------------------------------------------------------------------------
# clustering  : if true, simualted customers with the same preference patterns are merged on the fly
# seed        : seed provided to the random number generator
# max_samples : maximum number of realizations that can be placed in memory
#               at the same time during the sampling process (default 1e7)
##########################################################################################################
#a[ns/p,d] : 1 ⟺ simulated customer ns (clustering=false) or profile p (clustering=true) would be captured by location d
#ω[ns/p]   : weight of simulated customer ns (clustering=false) or profile p (clustering=true)
##########################################################################################################
function simulation_MNL(v, w, q, S; clustering=false, seed=0, max_samples=10000000)
    if(seed != 0) Random.seed!(seed) end
    (N,D) = size(v)
    
    #Divide the scenarios into groups to limit the memory usage (at most max_samples realizations of r.v. are simultaneously in memory)
    n_groups = Int(ceil(S*(D+1)*N/max_samples))
    if(S ≥ n_groups)
        group_sizes = ones(Int, n_groups) * Int(ceil(S/n_groups))
        group_sizes[end] = S-sum(group_sizes[1:end-1])
    else
        group_sizes = ones(Int, S)
    end
    
    #Initially, no observed preference patterns
    a=Bool[]
    ω=[]

    #Generate the scenarios by group and update the weight of the observed binary vectors
    for S_group in group_sizes
        #generate "S_group" scenarios
        ξ_precomputed = rand(Gumbel(0, 1), (N*S_group)*(D+1))

        #compute the binary vectors
        a, ω = compute_aω(v, w, q, ξ_precomputed, S_group, clustering=clustering, a_init=a, ω_init=ω)
    end
    
    #Sort the realizations of binary vectors in decreasing order of weight
    AΩ = hcat(a, ω)
    AΩ = AΩ[sortperm(AΩ[:, end], rev=true), :]
    a = Bool.(AΩ[:,1:end-1])
    ω = AΩ[:,end]
    
    #Return binary vectors a and their weights ω
    return (a, ω)
end

##########################################################################################################
######### Simulation for MMLN model in which a random parameter β[n]∼Unif[β_min, β_max], i.i.d. ##########
######################### multiplies the deterministic utilities of alternatives #########################
##########################################################################################################
# v[n,d] : deterministic utility of location d for customer n
# w[n,e] : deterministic utility of competitor e for customer n
# S      : number of scenarios
#---------------------------------------------------------------------------------------------------------
# seed        : seed provided to the random number generator
##########################################################################################################
# v_sim[r,d] : deterministic utility of location d for simulated customer r (S realizations per customer n)
# w_sim[r,e] : deterministic utility of competitor c for simulated customer r (S realizations per customer n)
##########################################################################################################
function simulation_β_MMNL(v, w, β_min, β_max; seed=0, S=1)
    if(seed != 0) Random.seed!(seed) end
    (N,D) = size(v)
    (N,E) = size(w)

    β = (β_max-β_min)*rand(N).+β_min
    v_sim = v.*β
    w_sim = w.*β
    for s ∈ 2:S
        β = (β_max-β_min)*rand(N).+β_min
        v_sim = vcat(v_sim, v.*β)
        w_sim = vcat(w_sim, w.*β)
    end

    return v_sim, w_sim
end

##########################################################################################################
############## Compute preference profiles a and their weights ω using precomputed utilities #############
##########################################################################################################
# v[n,d] : deterministic utility of potential location d for customer n
# w[n,c] : deterministic utility of competitor location c for customer n
# q[n]   : weight of customer n
# ξ      : precomputed vector of N*S×C realizations of errors terms 
# S      : number of scenarios
#---------------------------------------------------------------------------------------------------------
# clustering : if true, simualted customers with the same preference patterns are merged
# a_init     : previously computed binary vectors a
# ω_init     : previously computed weights ω
##########################################################################################################
# a : updated binary vectors a
# ω : updated weights ω
##########################################################################################################
function compute_aω(v, w, q, ξ, S; clustering=false, a_init=Bool[], ω_init=[])
    (N,D) = size(v)
    C = D+length(w[1,:])
    
    #Compute the behavior of the N*S simulated followers
    vr = repeat(v,S)
    wr = repeat(w,S)
    ω = repeat(q,S)
    ξr = reshape(ξ[1:N*S*C],(N*S,C))
    a = compute_a(vr, wr, ξr)
    
    #Add previous observations to the dataset
    if(a_init != [])
        a = vcat(a,a_init)
        ω = vcat(ω,ω_init)
    end
    
    if(clustering)
        #Return binary vectors a (a_p) and their weights ω (ω_p)
        return exact_clustering(a, ω)
    else
        #Return binary vectors a (a_ns) and their weights ω (q_ns)
        return (a, ω)
    end
end

##########################################################################################################
######################## Compute preference profiles a using precomputed utilities #######################
##########################################################################################################
#v[ns,d]   : deterministic utility of potential location d for customer n
#w[ns,c]   : deterministic utility of competitor location c for customer n
#ξ         : precomputed vector of N*S×C realizations of errors terms 
##########################################################################################################
#a[ns,d] : 1 ⟺ simulated customer ns is captured by location d
##########################################################################################################
function compute_a(v,w,ξ)
    (NS,D) = size(v)
    (NS,C) = size(ξ)
    u = v + ξ[:,1:D] #u_nd = total utility of alternative d for simulated customer ns
    W = maximum(w + ξ[:,D+1:C], dims=2) #W_ns = utility of the most attractive competing facility for ns
    a = Bool.(u.>W) 
    return a
end

##########################################################################################################
######################### Compute preference profiles a_p and their weights ω_p ##########################
###################### using precomputed simulated customers a_ns with weights q_ns ######################
##########################################################################################################
# a_sim[ns,d] : 1 ⟺ simulated customer ns is captured by location d
# q_sim[ns]   : weight of simulated customer ns
##########################################################################################################
# a[p,d] : 1 ⟺ preference profile p is captured by location d
# ω[p]   : weight of preference profile p
##########################################################################################################
function exact_clustering(a_sim, q_sim; remove_trivial=false)
    #Place the simulated customers sharing exactly the same preference profile in the same cluster
    D = size(a_sim)[2]
    rand_vals = rand(D)
    ID = a_sim*rand_vals
    IDreduced = unique(ID)
    a = a_sim[indexin(IDreduced,ID),:]
    idx_cluster = indexin(ID,IDreduced)
    n_clus = length(IDreduced)
    ω = zeros(n_clus)
    for i in 1:length(q_sim)
        ω[idx_cluster[i]] += q_sim[i]
    end

    #Remove the trivial pattern [0,...,0]
    if(remove_trivial)
        idx_non_trivial = findall(active->active≥1, sum(a,dims=2)[:,1])
        a = a[idx_non_trivial,:]
        ω = ω[idx_non_trivial]
    end
    
    #Sort the realizations of binary vectors in decreasing order of weight
    AΩ = hcat(a, ω)
    AΩ = AΩ[sortperm(AΩ[:, end], rev=true), :]
    a = Bool.(AΩ[:,1:end-1])
    ω = AΩ[:,end]
    
    #Return binary vectors a and their weights ω
    return (a, ω)
end