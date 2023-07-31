# Author: Robin Legault <legault@mit.edu>

#Packages
using Distributions

#Dependencies
include("Constants.jl")
include("IO_functions.jl")
include("Simulation_and_clustering_functions.jl")

##########################################################################################################
######################  Generate MNL instance with uniformly distributed customers, ###################### 
################################# as described in Haase and Müller (2014) ################################
##########################################################################################################
# N : number of customers
# D : number of available locations
# E : number of competing facilities
# instance : instance name e.g. "HM14_N_1000_D_100_E_10_α_1_β_1"
#---------------------------------------------------------------------------------------------------------
# width : width of the square on which customers are uniformly distributed
# seed        : seed provided to the random number generator
##########################################################################################################
function generate_instance_HM14(E, D, N, instance; width=30, seed=0)
    if(seed != 0) Random.seed!(seed) end
    coord_E = width*rand(E,2)
    coord_D = width*rand(D,2)
    coord_N = width*rand(N,2)
    c = zeros(N,D) #c[n,d] = L1 distance between customer n and location d
    k = zeros(N,E) #k[n,e] = L1 distance between customer n and competing facility e
    for n in 1:N
        for d in 1:D
            c[n,d] = sum(abs.(coord_N[n,:] - coord_D[d,:]))
        end
        for e in 1:E
            k[n,e] = sum(abs.(coord_N[n,:] - coord_E[e,:]))
        end
    end
    c = round.(c, digits = 3)
    k = round.(k, digits = 3)
    q  = Int.(ones(N))
    write_MNL_instance(instance, N, D, q, c, k; E=E)
end


##########################################################################################################
####################### Uniformly generate locations of different types in a square ######################
##########################################################################################################
# Dl[l] : number of available locations of type l
# El[l] : number of competing facilities of type l
#---------------------------------------------------------------------------------------------------------
# width : width of the square on which customers are uniformly distributed
# seed  : seed provided to the random number generator
##########################################################################################################
# coord_D    : coordinates of candidate locations
# coord_E    : coordinates of competing facilities
# type_D     : type_D[d]=l⟺ candidate location d is type l
# type_E     : type_E[e]=l ⟺ competing facility e is type l
##########################################################################################################
function generate_locations_mixture(Dl, El; width=30, seed=0)
    if(seed != 0) Random.seed!(seed) end

    D = sum(Dl)
    E = sum(El)
    L = length(Dl) #number of types of locations
    
    coord_D = width*(rand(D,2).-0.5)
    type_D = []
    for l in 1:L
        type_D = vcat(type_D, l*ones(Int,Dl[l]))
    end

    coord_E = width*(rand(E,2).-0.5)
    type_E = []
    for l in 1:L
        type_E = vcat(type_E, l*ones(Int,El[l]))
    end

    return coord_D, coord_E, type_D, type_E
end

##########################################################################################################
########## Generate customers of different types whose location follow a mixture of Gaussians ############
################### (this function uses fixed parameters, as presented in Section 6.2) ###################
##########################################################################################################
# N : number of customers
#---------------------------------------------------------------------------------------------------------
# seed  : seed provided to the random number generator
##########################################################################################################
# coord_N      : coordinates of customers
# type_N       : type[n]=t ⟺ customer n is type t
# q            : q[n] = weight of customer n
##########################################################################################################
function generate_population_mixture(N; seed=0)
    if(seed != 0) Random.seed!(seed) end

    # parameters of generative model (attributes Θ)
    # - three types of customers (e.g. teenagers, adults and seniors)
    # - four neighborhoods (e.g. downtown, working-class neighborhood, suburb 1, suburb 2)

    # π[i] = proportion of the population in neighborhood i
    π = [0.4, 0.3, 0.2, 0.1] #20% downtown, 30% working-class neighborhood, 40% suburb 1, 10% suburb 2

    # ρ[i,j] = proportion of customers of type j in neighborhood i
    ρ = [[0.2,0.7,0.1], #downtown                  : 20% teenagers, 70% adults, 10% seniors
         [0.3,0.4,0.3], #working-class neighborhood: 30% teenagers, 40% adults, 30% seniors
         [0.3,0.4,0.3], #suburb 1                  : 40% teenagers, 50% adults, 10% seniors
         [0.0,0.2,0.8]] #suburb 2                  :  0% teenagers, 20% adults, 80% seniors

    # spatial distribution of the customers in each neighborhood (bivariate Gaussian with parameters μ[i], Σ[i] for neighborhood i)
    μ = [[2.0,-2.0],   #downtown
         [-10,-10],    #working-class neighborhood
         [-4.0,10.0],  #suburb 1
         [12.0, -5.0]] #suburb 2

    Σ = [[[9.0 1.0];[1.0 9.0]],    #downtown
          [[9.0 -6.0];[-6.0 9.0]], #working-class neighborhood
          [[16.0 1.0];[1.0 4.0]],  #suburb 1
          [[2.0 0.0];[0.0 21.0]]]  #suburb 2
    
    # generate population
    q  = Int.(ones(N)) #same weight for each customer
    neighborhood = rand(Multinomial(1, π), N)'*[1,2,3,4] #neighborhood[n]=i ⟺ customer n is in neighborhood i
    type_N = zeros(Integer, N)
    for n=1:N
        type_N[n] = rand(Multinomial(1, ρ[neighborhood[n]]))'*[1,2,3] #type_N[n]=t ⟺ customer n is type t
    end
    coord_N = zeros(N,2) #coordinates of customers
    for n=1:N
        coord_N[n,:] = rand(MvNormal(μ[neighborhood[n]], Σ[neighborhood[n]]))
    end
    
    return coord_N, type_N, q
end

##########################################################################################################
################ Generate the utilities for a set of customers N and locations C=D∪E #####################
########################## according to the MMNL model presented in Section 6.2 ##########################
################### (this function uses fixed parameters, as presented in Section 6.2) ###################
##########################################################################################################
# α       : aversion to competing facilities
# coord_D : coordinates of candidate locations
# coord_E : coordinates of competing facilities
# type_D  : type_D[d]=l⟺ candidate location d is type l
# type_E  : type_E[e]=l ⟺ competing facility e is type l
# coord_N : coordinates of customers
# type_N  : type[n]=t ⟺ customer n is type t
#---------------------------------------------------------------------------------------------------------
# seed  : seed provided to the random number generator
##########################################################################################################
# c[n,d] : distance from customer n to location d
# k[n,e] : distance from customer n competing facility e
# v[n,d] : deterministic utility of location d for customer n
# w[n,e] : deterministic utility of competing facility e for customer n
##########################################################################################################
function generate_utilities_mixture(α, coord_D, coord_E, type_D, type_E, coord_N, type_N; seed=0)
    if(seed != 0) Random.seed!(seed) end
    
    N = length(type_N)
    D = length(type_D)
    E = length(type_E)
    
    # γ[t,l] = aversion of customers of type t for locations of type l
    γ = [[20,60,40], #teenagers : (separate shop) prefered to (mall) prefered to (rest area) 
         [40,20,60], #adults    : (rest area) prefered to (separate shop) prefered to (mall)
         [60,40,20]] #seniors   : (mall) prefered to (rest area) prefered to (separate shop)

    # τ[t] = aversion of customers of type t to distances
    τ = [3,1,2] # teenagers are more distance-averse than seniors than adults
    
    # compute distances
    v = zeros(N,D)
    w = zeros(N,E)
    c = zeros(N,D)
    k = zeros(N,E)
    for n in 1:N
        for j in 1:D
            c[n,j] = sum(abs.(coord_N[n,:] - coord_D[j,:]))
        end
        for l in 1:E
            k[n,l] = sum(abs.(coord_N[n,:] - coord_E[l,:]))
        end
    end

    # utility of candidate locations
    for n=1:N
        for d=1:D
            v[n,d] -= τ[type_N[n]]*c[n,d]     # disutility associated with travel distance
            v[n,d] -= γ[type_N[n]][type_D[d]] # disutility associated with type of location
        end
    end

    # utility of competing facilities
    for n=1:N
        for e=1:E
            w[n,e] -= τ[type_N[n]]*k[n,e]     # disutility associated with travel distance
            w[n,e] -= γ[type_N[n]][type_E[e]] # disutility associated with type of location
        end
    end
    w = α.*w # competitor aversion multiplier
    v = round.(v, digits = 3)
    w = round.(w, digits = 3)

    return v, w
end

##########################################################################################################
########### Generate an instance for a set of customers N and locations C=D∪E according to the ###########
####  model presented in Section 6.2. The realizations of β in the MMNL model are not generated here. ####
################### (this function uses fixed parameters, as presented in Section 6.2) ###################
##########################################################################################################
# N        : number of customers
# coord_D : coordinates of candidate locations
# coord_E : coordinates of competing facilities
# type_D  : type_D[d]=l⟺ candidate location d is type l
# type_E  : type_E[e]=l ⟺ competing facility e is type l
# α        : aversion to competing facilities
# instance : instance name e.g. "MIX_N_10000_D1_15_D2_15_D3_15_E1_5_E2_5_E3_5_α_1"
# β_min    : minimum multiplier β (β[n]∼Unif[β_min, β_max], i.i.d.)
# β_max    : maximum multiplier β
#---------------------------------------------------------------------------------------------------------
# seed : seed provided to the random number generator
# S    : number of realizations of the random coefficients β in the MMNL model
##########################################################################################################
function generate_instance_MMNL_mixture(N, coord_D, coord_E, type_D, type_E, α, instance, β_min, β_max; seed=0, S=1)
    if(seed != 0) Random.seed!(seed) end

    D = length(type_D)
    E = length(type_E)
    coord_N, type_N, q = generate_population_mixture(N)
    v, w = generate_utilities_mixture(α, coord_D, coord_E, type_D, type_E, coord_N, type_N)
    v, w = simulation_β_MMNL(v, w, β_min, β_max; S=S)

    # Save attributes Θ of the customers
    write_MNL_instance_utilities(instance, N, D, q, v, w; E=E)
end