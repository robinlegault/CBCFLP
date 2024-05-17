# Author: Robin Legault <legault@mit.edu>

#Packages
using LinearAlgebra

##########################################################################################################
############################ Objective value of solution x under the MNL model ###########################
##########################################################################################################
# x[d]       : 1 ⟺ location d is open in the solution
# q[n]       : weight of customer n
# V[n,d]     : exp(v[n,d]), where v[n,d] is the deterministic utility of location d for customer n
# W_super[n] : exp(w_super[n]), where w_super[n] is the deterministic utility of the competing facility for customer n
##########################################################################################################
# Z_MNL : Objective value of solution x under the MNL model
##########################################################################################################
function objective_value_MNL(x,q,V,W_super)
    Z_MNL = sum(q)-dot(q,W_super./(V*x+W_super))
    return Z_MNL
end

##########################################################################################################
############### Probability for each customer to select each alternative under the MNL model #############
##########################################################################################################
# x[d]   : 1 ⟺ location d is open in the solution
# v[n,d] : exp(v[n,d]), where v[n,d] is the deterministic utility of location d for customer n
# w[n]   : exp(w_super[n]), where w_super[n] is the deterministic utility of the competing facility for customer n
##########################################################################################################
# prob_MNL[n,c] : Probability for customer n to select location c ∈ D∪E
##########################################################################################################
function probabilities_MNL(x,v,w)
    (N,D) = size(v)
    vx = v.-(round.(ones(D)-x)*1e10)' #infinite disutility for closed alternatives
    vw = hcat(vx,w)
    vw .-= maximum(vw,dims=2)
    prob_MNL = exp.(vw) ./ sum(exp.(vw), dims=2)
    return prob_MNL
end