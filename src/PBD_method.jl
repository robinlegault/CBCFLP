# Author: Robin Legault <legault@mit.edu>

#Packages
using CPLEX, JuMP

#Dependencies
include("Greedy_method.jl")
include("Entropy_functions.jl")
include("Data_structures.jl")

##########################################################################################################
################################# PBD method for solving model \hat{DEQ} #################################
##########################################################################################################
# a[p,d]: 1 ⟺ profile p is captured by location d
# ω[p]  : weight of profile p
# b     : budget
#---------------------------------------------------------------------------------------------------------
# min_P_submod  : If |P|<min_P_submod, set P1=P. Using submodular cuts is generally not advantageous for very small instances
# tvd_min       : If δ*<tvd_min, set P1=P. Using submodular cuts is generally not advantageous when the total variation distance is small
# P1_auto       : if true, P1 is selected using the knee method (PBD)
# P1_weight_min : if P1_auto=false, P1 is smallest set s.t. Ω ≥ P1_weight_min. 0->Pure submodular, 1->Pure 0-1
# heuristic     : if true, a greedy solution is first computed and provided to the solver
# display       : if 0, no output is printed during the solving process. Higher values provide more details
# limit_time    : time limit provided to the solver, in seconds
# stats_PBD     : mutable structure Stats_PBD initially containing details on the instance
# return_stats  : if true, return populated Stats_PBD structure. Otherwise, return opt. sol and value
##########################################################################################################
# stats_PBD: updated with the solution and computational statistics (returned if return_stats)
# x_hat[d] : 1 ⟺ location d is open; optimal solution of model \hat{DEQ} (returned if !return_stats)
# Z_hat    : optimal value of model \hat{DEQ}                             (returned if !return_stats)
##########################################################################################################
function PBD(a, ω, b; min_P_submod=0, tvd_min=0, P1_auto=true, P1_weight_min=1, heuristic=true, display=0, limit_time=3600, stats_PBD=[], return_stats=true)  

    #statistics
    (P,D) = size(a)
    if(stats_PBD==[])
        stats_PBD = init_Stats_PBD(D)
    end

    #If |P|<min_P_submod, set P1=P
    if(P < min_P_submod)
        P1_auto = false
        P1_weight_min=1
    end

    #Sort the realizations of binary vectors in decreasing order of weight, if necessary
    if(!issorted(reverse(ω)))
        AΩ = hcat(a, ω)
        AΩ = AΩ[sortperm(AΩ[:, end], rev=true), :]
        a = Bool.(AΩ[:,1:end-1])
        ω = AΩ[:,end]
    end
    
    #Build P1
    cum_weights = cumsum(ω)./sum(ω)
    if(P1_auto) #P1 is selected using the knee method
        stats_PBD.P1_auto = true
        x = collect(1:P)./P
        y = cum_weights
        diffs = y-x
        P1 = findlast(i -> i ==  maximum(diffs), diffs)
        tvd = diffs[P1]
        stats_PBD.δ = tvd
        #  If δ*<tvd_min, set P1=P
        if(tvd < tvd_min)
            P1 = P
        end
    else #Relative weight of P1 is given, P1 is the smallest set of profiles that represent at least this weight
        stats_PBD.P1_weight_min = P1_weight_min
        if(P1_weight_min==0)
            P1 = 0
        elseif(P1_weight_min==1)
            P1 = P
        else
            P1 = findfirst(weight -> weight ≥ P1_weight_min, cum_weights)
        end
    end
    if(P1>0)
        stats_PBD.P1_weight = cumsum(ω)[P1]/sum(ω)
    else
        stats_PBD.P1_weight = 0.0
    end
    P2 = P-P1
    stats_PBD.P = P
    stats_PBD.P1 = P1
    stats_PBD.P2 = P2
    
    #Preference profiles and weights for sets P1 and P2=P\P1
    a1 = a[1:P1,:]
    ω1 = ω[1:P1]
    a2 = a[P1+1:end,:]
    ω2 = ω[P1+1:end]

    #Model
    m = Model()

    #CPLEX parameters
    if(display==0) set_silent(m) end
    set_optimizer(m, CPLEX.Optimizer)
    if(display≥1) set_optimizer_attribute(m, "CPX_PARAM_MIPDISPLAY", display) end # Display model details
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", limit_time) # Time limit

    #Variables
    @variable(m, x[d=1:D], Bin) #xd=1 ⟺ location d is open
    @variable(m, y[p=1:P1] ≤ 1) #yp=1 ⟺ the demand of profile p is captured
    @variable(m, Θ ≤ sum(ω2))   #initial UB on the captured weight on set P2: weight of set P2

    #Facility location constraints
    @constraint(m, con_x, sum(x[d] for d in 1:D) == b) #at most b location can be opened (∃ opt. sol. with exactly b open locations)

    #Demand capture constraints
    @constraint(m, con_y[p=1:P1], y[p] ≤ sum(a1[p,d]*x[d] for d=1:D)) #y[p]=0 if no location d∈D(p) is open

    #Objective: Maximize the weighted captured demand
    @objective(m, Max, sum(ω1[p]*y[p] for p=1:P1) + Θ)
    
    #Greedy initial solution
    if(heuristic)
        #Greedy heuristic
        stats_PBD.time_greedy = (@timed begin (stats_PBD.x_greedy, stats_PBD.Z_greedy) = greedy_add(a, ω, b) end).time

        #Greedy solution captures all the demand -> optimal
        if(stats_PBD.Z_greedy == sum(ω))
            stats_PBD.x_opt = stats_PBD.x_greedy
            stats_PBD.Z_opt = sum(ω)
            stats_PBD.is_optimal = true
            return stats_PBD
        end

        #Feed initial solution to the model
        x_init = stats_PBD.x_greedy
        y_init = (a1 * x_init) .≥ 1
        Θ_init = ((a2 * x_init).≥1)'*ω2
        for d in 1:D
            set_start_value(x[d], x_init[d])
        end
        for p in 1:P1
            set_start_value(y[p], y_init[p])
        end
        set_start_value(Θ, Θ_init)
    end
    
    #Callback function: submodularity cuts (lazy constraints, called on integer nodes)
    if(P2 > 0)
        aω2 = a2.*ω2
        margin_gain_all_open = sum(max.(0, 2*aω2 .- sum(aω2, dims=2)), dims=1)[1,:]
        margin_gain_all_closed = sum(aω2,dims=1)[1,:]
        function callback_submodular_cuts(cb_data)
            stats_PBD.n_integer_nodes_BnC =  stats_PBD.n_integer_nodes_BnC + 1
            x_val = [Bool(round(callback_value(cb_data, x[d]))) for d in 1:D]
            Θ_val = callback_value(cb_data, Θ)
            active = findall(x->x==1,x_val)
            inactive = setdiff(collect(1:D), active)
            capture_P2 = sum(min.(aω2*x_val,ω2))
            if(Θ_val > capture_P2)
                stats_PBD.n_submodular_cuts = stats_PBD.n_submodular_cuts + 2
                margin_gain = sum(max.(0, aω2[:,inactive].-sum(aω2[:,active], dims=2)), dims=1)[1,:]
                margin_loss = sum(min.(a2[:,active], a2[:,active].==sum(a2[:,active], dims=2)[:,1]).*ω2, dims=1)[1,:]
                cut1 = @build_constraint(Θ <= capture_P2 + sum(margin_gain[i]*x[d] for (i,d) in enumerate(inactive)) - sum(margin_gain_all_open[d]*(1-x[d]) for d in active)) #(27)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cut1)
                cut2 = @build_constraint(Θ <= capture_P2 + sum(margin_gain_all_closed[d]*x[d] for d in inactive) - sum(margin_loss[i]*(1-x[d]) for (i,d) in enumerate(active))) #(37)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cut2)
            end
        end
        set_attribute(m, MOI.LazyConstraintCallback(), callback_submodular_cuts)
    end
    
    #Solve model \hat{DEQ}
    optimize!(m)
    stats_PBD.time_exact = solve_time(m)

    # Get the number of visited nodes
    stats_PBD.n_nodes_BnC = MOI.get(m, MOI.NodeCount())
    
    #Return the optimal solution, value and other statistics in struct stats_PBD
    x_opt = Bool.(round.(value.(x)))
    Z_opt = objective_value(m)
    stats_PBD.x_opt = x_opt
    stats_PBD.Z_opt = Z_opt
    stats_PBD.is_optimal = (termination_status(m) == MOI.OPTIMAL)

    if(return_stats)
        return stats_PBD
    else
        return (x_opt, Z_opt)
    end
end