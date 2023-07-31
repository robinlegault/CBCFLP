# Author: Robin Legault <legault@mit.edu>

##########################################################################################################
################################## Greedy heuristic for model \hat{DEQ} ##################################
##########################################################################################################
# a[p,d] : 1 ⟺ profile p is captured by location d
# ω[p]   : weight of profile p
# b      : budget
##########################################################################################################
# x_greedy[d] : 1 ⟺ location d is open; greedy solution for model \hat{DEQ}
# Z_greedy    : LB on the optimal value of model \hat{DEQ}
##########################################################################################################
function greedy_add(a, ω, b)
    D = size(a)[2]
    aω = a.*ω
    obj = 0
    rem_idx = collect(1:D)
    for iter in 1:b
        gains = sum(aω, dims=1)[1,:]
        obj += maximum(gains) 
        d_star = argmax(gains) #add the remaining facility with the largest gain in capture
        deleteat!(rem_idx, d_star)
        aω = aω[findall(i -> aω[i,d_star]==0, collect(1:size(aω)[1])),setdiff(1:end, d_star)] #remove the captured profiles
    end
    x_greedy = zeros(Bool, D)
    x_greedy[setdiff(1:D, rem_idx)].=1
    Z_greedy = obj

    return (x_greedy, Z_greedy)
end
