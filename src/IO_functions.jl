# Author: Robin Legault <legault@mit.edu>

#Packages
using XLSX, DelimitedFiles

#Dependencies
include("Constants.jl")

##########################################################################################################
##### Read MNL conditional instances files (.data), parameters as described in Mai and Lodi (2020)  ######
##########################################################################################################
# instance : instance name e.g. "NYC_82341_59_5_1"
##########################################################################################################
# N      : number of customers
# D      : number of available facility locations
# E      : number of competing facilities
# q[n]   : weight of customer n
# c[n,d] : distance from customer n to location d
# k[n,e] : distance from customer n competing facility e
##########################################################################################################
function read_MNL_instance(instance)
    file_name = INSTANCES_PATH*instance*".data"
    inst_file = open(file_name) do file read(file, String) end
    inst_file = split(inst_file, "\n")
    line = parse.(Int, split(inst_file[1], " "))
    N = line[1] 
    D = line[2]
    E = line[3] 
    q = parse.(Float64, split(inst_file[2], " ")) 
    k = zeros(N,E)
    for i in 3:2+E
        k[:,i-2] = -parse.(Float64, split(inst_file[i], " "))
    end
    c = zeros(N,D)
    for i in 3+E:2+E+D
        c[:,i-(2+E)] = -parse.(Float64, split(inst_file[i], " "))
    end
    return N, D, E, q, c, k
end

##########################################################################################################
############# Read MNL conditional instances files (.data) storing deterministic utilities  ##############
##########################################################################################################
# instance : instance name e.g. "NYC_82341_59_5_1"
##########################################################################################################
# N          : number of customers
# D          : number of available facility locations
# E          : number of competing facilities
# q[n]       : weight of customer n
# v[n,d]     : deterministic utility of location d for customer n
# V[n,d]     : exp(v[n,d])
# w_super[n] : deterministic utility of the unique (super-)competing alternative for customer n
# W_super[n] : exp(w_super[n])
##########################################################################################################
function read_MNL_instance_utilities(instance)
    file_name = INSTANCES_PATH*"utilities_"*instance*".data"
    inst_file = open(file_name) do file read(file, String) end
    inst_file = split(inst_file, "\n")
    line = parse.(Int, split(inst_file[1], " "))
    N = line[1] 
    D = line[2]
    E = line[3] 
    q = parse.(Float64, split(inst_file[2], " ")) 
    w = zeros(N,E)
    for i in 3:2+E
        w[:,i-2] = parse.(Float64, split(inst_file[i], " "))
    end
    v = zeros(N,D)
    for i in 3+E:2+E+D
        v[:,i-(2+E)] = parse.(Float64, split(inst_file[i], " "))
    end

    W_super = sum(exp.(w),dims=2)
    w_super = log.(W_super)

    #rescale utilities for numerical stability (does not modify the instance)
    v = v.-w_super      
    w_super = zeros(N) 

    V  = exp.(v)
    W_super = exp.(w_super)
    return N, D, E, q, v, w_super, V, W_super
end

##########################################################################################################
##### Write MNL conditional instances files (.data), parameters as described in Mai and Lodi (2020)  #####
##########################################################################################################
# instance : instance name e.g. "NYC_82341_59_5_1"
# N        : number of customers
# D        : number of available facility locations
# q[n]     : weight of customer n
# c[n,d]   : distance from customer n to location d
# k[n,e] : distance from customer n competing facility e
#---------------------------------------------------------------------------------------------------------
# E : number of competing facilities
# folder : write the instance to a subfolder inside INSTANCES_PATH, e.g. "entropy/"
##########################################################################################################
function write_MNL_instance(instance, N, D, q, c, k; E=1, folder="")
    file_name = INSTANCES_PATH*folder*instance*".data"
    open(file_name, "w") do io
        writedlm(io, [N;D;E]', " ")
        writedlm(io, Integer.(q)', " ")
        for e in 1:E
            writedlm(io, -k[:,e]', " ")
        end
        for d in 1:D
            writedlm(io, -c[:,d]', " ")
        end
    end
end

##########################################################################################################
########### Write MNL conditional instances files (.data) by storing deterministic utilities  ############
##########################################################################################################
# instance : instance name e.g. "NYC_82341_59_5_1"
# N        : number of customers
# D        : number of available facility locations
# q[n]     : weight of customer n
# v[n,d]   : deterministic utility of location d for customer n
# w[n,e]   : deterministic utility of the competing facility e for customer n
#---------------------------------------------------------------------------------------------------------
# E : number of competing facilities
# folder : write the instance to a subfolder inside INSTANCES_PATH, e.g. "entropy/"
##########################################################################################################
function write_MNL_instance_utilities(instance, N, D, q, v, w; E=1, folder="")
    file_name = INSTANCES_PATH*folder*"utilities_"*instance*".data"
    open(file_name, "w") do io
        writedlm(io, [N;D;E]', " ")
        writedlm(io, Integer.(q)', " ")
        for e in 1:E
            writedlm(io, w[:,e]', " ")
        end
        for d in 1:D
            writedlm(io, v[:,d]', " ")
        end
    end
end

##########################################################################################################
###### Parse MNL conditional instances from file, parameters as described in Mai and Lodi (2020)  ########
##########################################################################################################
# instance : instance name e.g. "NYC_82341_59_5_1"
# α        : aversion to competing facilities
# β        : aversion to distances
#---------------------------------------------------------------------------------------------------------
# α0      : aversion to competing facilities (additive penalty)
##########################################################################################################
# N          : number of customers
# D          : number of available facility locations
# E          : number of competing facilities
# q[n]       : weight of customer n
# v[n,d]     : deterministic utility of location d for customer n
# V[n,d]     : exp(v[n,d])
# w_super[n] : deterministic utility of the unique (super-)competing alternative for customer n
# W_super[n] : exp(w_super[n])
##########################################################################################################
function parse_MNL_instance(instance, α, β; α0=0)
    N, D, E, q, c, k = read_MNL_instance(instance)
    v  = -β.*c 
    w  = -(α*β.*k).-α0
    W_super = sum(exp.(w),dims=2)
    w_super = log.(W_super)

    #rescale utilities for numerical stability (does not modify the instance)
    v = v.-w_super      
    w_super = zeros(N) 

    V  = exp.(v)
    W_super = exp.(w_super)
    return N, D, E, q, v, w_super, V, W_super
end

##########################################################################################################
################## Parse conditional MNL instance from a set of location and parameters ##################
##########################################################################################################
# coord_D : coordinates of candidate locations
# coord_E : coordinates of competing facilities
# coord_N : coordinates of customers
# α       : aversion to competing facilities
# β       : aversion to distances
#---------------------------------------------------------------------------------------------------------
# α0      : aversion to competing facilities (additive penalty)
##########################################################################################################
# N          : number of customers
# D          : number of available facility locations
# E          : number of competing facilities
# q[n]       : weight of customer n
# v[n,d]     : deterministic utility of location d for customer n
# V[n,d]     : exp(v[n,d])
# w_super[n] : deterministic utility of the unique (super-)competing alternative for customer n
# W_super[n] : exp(w_super[n])
##########################################################################################################
function parse_MNL_instance_coord(coord_D, coord_E, coord_N, α, β; α0=0)
    D = size(coord_D)[1]
    N = size(coord_N)[1]
    E = size(coord_E)[1]
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
    q  = ones(N)
    v  = -β.*c
    w  = -(α*β.*k).-α0
    W_super = sum(exp.(w),dims=2)
    w_super = log.(W_super)
    
    #rescale utilities for numerical stability (does not modify the instance)
    v = v.-w_super      
    w_super = zeros(N) 

    V  = exp.(v)
    W_super = exp.(w_super)

    return N, D, E, q, v, w_super, V, W_super
end

##########################################################################################################
##################################### Append a row to an Excel file  #####################################
##########################################################################################################
# workbook_path  : outfile name of the excel file, e.g. "knee.xlsx"
# row_data       : data to append
#---------------------------------------------------------------------------------------------------------
# sheet_name     : append the row at the end of this sheet
# column_names   : name of the columns of the data 
# overwrite_last : if true, the appended row overwrites the last row in the file
##########################################################################################################
function append_xl_row(outfile, row_data; sheet_name="Sheet1", column_names=[], overwrite_last=false)
    
    workbook_path = RESULTS_OUTPUT_PATH*outfile #location: ../results folder
    if !isfile(workbook_path)
        # File doesn't exist, so create it using ExcelFiles.jl
        if(column_names==[])
            for i in 1:length(row_data)
                push!(column_names, "col$i")
            end
        end
        XLSX.writetable(workbook_path, zeros(length(column_names)), column_names; sheetname=sheet_name)
        overwrite_last = true
    end
    XLSX.openxlsx(workbook_path, mode="rw") do xf
        sheet = xf[sheet_name]
        num_rows = XLSX.get_dimension(sheet).stop.row_number
        
        if(overwrite_last)
            sheet[num_rows,:] = row_data
        else
            sheet[num_rows+1,:] = row_data
        end
    end
end

##########################################################################################################
#################### Append a row containing the fields of Stats_PBD to an Excel file ####################
##########################################################################################################
# stats_PBD  : populated structure Stats_PBD 
# row_data   : data to append
#---------------------------------------------------------------------------------------------------------
# sheet_name     : append the row at the end of this sheet
# column_names   : name of the columns of the data 
# overwrite_last : if true, the appended row overwrites the last row in the file
##########################################################################################################
function save_results_PBD(stats_PBD, outfile; sheet_name="Sheet1")
    row_stats = [            
            stats_PBD.inst_name               
            stats_PBD.D
            stats_PBD.b
            stats_PBD.N
            stats_PBD.S
            stats_PBD.H
            stats_PBD.Φ
            stats_PBD.δ
            stats_PBD.weight_trivial
            stats_PBD.weight_nontrivial
            stats_PBD.weight_total
            stats_PBD.time_clustering
            stats_PBD.time_greedy
            stats_PBD.time_exact
            stats_PBD.n_nodes_BnC
            stats_PBD.n_integer_nodes_BnC
            stats_PBD.n_submodular_cuts
            stats_PBD.P
            stats_PBD.P1
            stats_PBD.P2
            stats_PBD.P1_auto
            stats_PBD.P1_weight_min
            stats_PBD.P1_weight
            string(stats_PBD.x_greedy)         
            stats_PBD.Z_greedy
            stats_PBD.Z_greedy_stoch
            string(stats_PBD.x_opt)                   
            stats_PBD.Z_opt
            stats_PBD.Z_opt_stoch
            stats_PBD.limit_time
            stats_PBD.is_optimal
    ]
    append_xl_row(outfile, row_stats; column_names=column_names_PBD, sheet_name=sheet_name)
end

##########################################################################################################
######################### Print important fields of Stats_PBD to standard output #########################
##########################################################################################################
# stats_PBD  : populated structure Stats_PBD 
##########################################################################################################
function print_results_PBD(stats_PBD)
    if(stats_PBD.P1_auto)
        println("Method            : PBD, Ω = Ω_i*")
    elseif(stats_PBD.time_clustering == 0)
        println("Method            : SAA")
    else
        println("Method            : PBD, Ω = ",stats_PBD.P1_weight_min)
    end
    println("Instance          : ",stats_PBD.inst_name)
    println("Optimal           : ",stats_PBD.is_optimal)
    println("Weight total      : ",stats_PBD.weight_total/stats_PBD.S)
    println("Weight nontrivial : ",stats_PBD.weight_nontrivial/stats_PBD.S)
    println("Value             : ",stats_PBD.Z_opt,
    " (",round(100*stats_PBD.S*stats_PBD.Z_opt/stats_PBD.weight_nontrivial, digits=1),"% of nontrivial), (",
         round(100*stats_PBD.S*stats_PBD.Z_opt/stats_PBD.weight_total, digits=1),"% of total)")
    println("Value (stoch)     : ",stats_PBD.Z_opt_stoch)
    println("H                 : ",stats_PBD.H)
    println("δ*                : ",stats_PBD.δ)
    println("P1                : ",stats_PBD.P1)
    println("P2                : ",stats_PBD.P2)
    println("n_nodes_BnC       : ",stats_PBD.n_nodes_BnC)
    println("time_exact        : ",stats_PBD.time_exact,"\n")
end