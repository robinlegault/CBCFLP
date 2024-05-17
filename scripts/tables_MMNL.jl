# Author: Robin Legault <legault@mit.edu>

# Packages
using DataFrames, XLSX, Statistics, Plots, LaTeXStrings, Random

# Dependencies
include("../src/IO_functions.jl")
include("../src/MNL_functions.jl")

##########################################################################################################
#################### Computational results for generative MMNL instances, Section 6.2. ###################
##########################################################################################################

# Load results
Dl = [25, 25, 25]
El = [10, 10, 10]
α = 1
results_file = "MMNL_MIX_paper.xlsx" #
df = DataFrame(XLSX.readtable(RESULTS_PATH*results_file, "Sheet1"; header=true, infer_eltypes=true))
df_MOA = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*"GGX_MOA_"*results_file, "MOA"; header=true, infer_eltypes=true))
df_GGX = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*"GGX_MOA_"*results_file, "GGX"; header=true, infer_eltypes=true))
col_βmin_GGX = parse.(Float64, [match(r"βmin_([\d.]+)", string).captures[1] for string in df_GGX.inst_name])
col_βmin_MOA = parse.(Float64, [match(r"βmin_([\d.]+)", string).captures[1] for string in df_MOA.inst_name])
col_βmin = parse.(Float64, [match(r"βmin_([\d.]+)", string).captures[1] for string in df.inst_name])
col_βmax = parse.(Float64, [match(r"βmax_([\d.]+)", string).captures[1] for string in df.inst_name])
col_α = parse.(Float64, [match(r"α_([\d.]+)", string).captures[1] for string in df.inst_name])
df_GGX[!, "βmin"] = col_βmin_GGX
df_MOA[!, "βmin"] = col_βmin_MOA
df[!, "βmin"] = col_βmin
df[!, "βmax"] = col_βmax
df[!, "α"] = col_α
list_N = unique(df.N)
list_N_GGX = unique(df_GGX.N)
list_N_MOA = unique(df_MOA.N)
list_βmin = reverse(sort!(unique(df.βmin)))
list_βmax = reverse(sort!(unique(df.βmax)))
list_S = unique(df.S)
list_b = unique(df.b)
nβ = length(list_βmin)
nS = length(list_S)
nb = length(list_b)
nN = length(list_N)
nN_GGX = length(list_N_GGX)
nN_MOA = length(list_N_MOA)
global N_reference = 1000000

# Collect statistics
time_exact_PBD         = zeros(nβ, nN, nS)
time_exact_SAA         = zeros(nβ, nN, nS)
time_exact_01          = zeros(nβ, nN, nS)
time_total_PBD         = zeros(nβ, nN, nS)
time_total_SAA         = zeros(nβ, nN, nS)
time_total_01          = zeros(nβ, nN, nS)
time_MOA               = zeros(nβ, nN_MOA)
time_GGX               = zeros(nβ, nN_GGX)
tlim_reached_PBD       = zeros(nβ, nN, nS)
tlim_reached_SAA       = zeros(nβ, nN, nS)
tlim_reached_01        = zeros(nβ, nN, nS)
list_H_hat             = zeros(nβ, nN, nS)
list_H_hat_avg         = zeros(nβ, nS)
list_n_submodular_cuts = zeros(nβ, nN, nS)
list_n_nodes_BnC_PBD   = zeros(nβ, nN, nS)
list_n_nodes_BnC_SAA   = zeros(nβ, nN, nS)
list_n_nodes_BnC_01    = zeros(nβ, nN, nS)
list_P                 = zeros(nβ, nN, nS)
list_P1                = zeros(nβ, nN, nS)
list_δ                 = zeros(nβ, nN, nS)
list_RGap_MOA          = zeros(nβ, nN_MOA)
list_RGap_GGX          = zeros(nβ, nN_GGX)
list_RGap_SIM          = zeros(nβ, nN, nS)
list_Z_SIM             = zeros(nβ, nN, nS)
list_Z_MNL_MOA         = zeros(nβ, nN_MOA)
list_Z_MNL_GGX         = zeros(nβ, nN_GGX)
list_Z_MNL_SIM         = zeros(nβ, nN, nS)
list_Z_MMNL_MOA        = zeros(nβ, nN_MOA)
list_Z_MMNL_GGX        = zeros(nβ, nN_GGX)
list_Z_MMNL_SIM        = zeros(nβ, nN, nS)
list_RGap_MOA_b        = zeros(nβ, nN_MOA, nb)
list_RGap_GGX_b        = zeros(nβ, nN_GGX, nb)
list_RGap_SIM_b        = zeros(nβ, nN, nS, nb)
list_Z_SIM_b           = zeros(nβ, nN, nS, nb)
list_Z_MNL_MOA_b       = zeros(nβ, nN_MOA, nb)
list_Z_MNL_GGX_b       = zeros(nβ, nN_GGX, nb)
list_Z_MNL_SIM_b       = zeros(nβ, nN, nS, nb)
list_Z_MMNL_MOA_b      = zeros(nβ, nN_MOA, nb)
list_Z_MMNL_GGX_b      = zeros(nβ, nN_GGX, nb)
list_Z_MMNL_SIM_b      = zeros(nβ, nN, nS, nb)
list_params            = Array{String}(undef, nβ, nN, nS)
list_params_GGX        = Array{String}(undef, nβ, nN_GGX)
list_params_MOA        = Array{String}(undef, nβ, nN_MOA)
for iβ ∈ 1:nβ
    βmin = list_βmin[iβ]
    βmax = list_βmax[iβ]
    instance = "MIX_N_"*string(N_reference)*"_D1_"*string(Dl[1])*"_D2_"*string(Dl[2])*"_D3_"*string(Dl[3])*"_E1_"*string(El[1])*
                    "_E2_"*string(El[2])*"_E3_"*string(Dl[3])*"_α_"*string(α)*"_βmin_"*string(βmin)*"_βmax_"*string(βmax)
    trash, D, E, q, v, w_super, V, W_super = read_MNL_instance_utilities(instance)
    df_β = filter(row -> row.βmin==βmin, df)
    df_MOA_β = filter(row -> row.βmin==βmin, df_MOA)
    df_GGX_β = filter(row -> row.βmin==βmin, df_GGX)
    for iN ∈ 1:nN
        N = list_N[iN]
        df_βN = filter(row -> row.N==N, df_β)
        df_PBD_βN = filter(row -> (row.P1_auto), df_βN)
        df_SAA_βN  = filter(row -> (row.time_clustering==0.0), df_βN)
        df_01_βN  = filter(row -> (row.time_clustering>0.0 && !row.P1_auto && row.P1_weight_min==1), df_βN)
        df_PBD_βN_opt = filter(row -> row.is_optimal, df_PBD_βN)
        df_SAA_βN_opt  = filter(row -> row.is_optimal, df_SAA_βN)
        df_01_βN_opt  = filter(row -> row.is_optimal, df_01_βN)

        for iS ∈ 1:nS
            S = list_S[iS]
            df_βNS = filter(row -> (row.S==S), df_βN)
            df_PBD_βNS = filter(row -> (row.S==S), df_PBD_βN)
            df_SAA_βNS  = filter(row -> (row.S==S), df_SAA_βN)
            df_01_βNS  = filter(row -> (row.S==S), df_01_βN)
            df_PBD_βNS_opt = filter(row -> (row.S==S), df_PBD_βN_opt)
            df_SAA_βNS_opt  = filter(row -> (row.S==S), df_SAA_βN_opt)
            df_01_βNS_opt  = filter(row -> (row.S==S), df_01_βN_opt)

            list_n_nodes_BnC_PBD[iβ,iN,iS]   = mean(df_PBD_βNS_opt.n_nodes_BnC)
            list_n_nodes_BnC_SAA[iβ,iN,iS]   = mean(df_SAA_βNS_opt.n_nodes_BnC)
            list_n_nodes_BnC_01[iβ,iN,iS]    = mean(df_01_βNS_opt.n_nodes_BnC)
            time_exact_PBD[iβ,iN,iS] = mean(df_PBD_βNS_opt.time_exact)
            time_exact_SAA[iβ,iN,iS] = mean(df_SAA_βNS_opt.time_exact)
            time_exact_01[iβ,iN,iS]  = mean(df_01_βNS_opt.time_exact)
            time_total_PBD[iβ,iN,iS] = mean(df_PBD_βNS_opt.time_clustering) + mean(df_PBD_βNS_opt.time_greedy) + mean(df_PBD_βNS_opt.time_exact)
            time_total_SAA[iβ,iN,iS] = mean(df_SAA_βNS_opt.time_clustering) + mean(df_SAA_βNS_opt.time_greedy) + mean(df_SAA_βNS_opt.time_exact)
            time_total_01[iβ,iN,iS]  = mean(df_01_βNS_opt.time_clustering) + mean(df_01_βNS_opt.time_greedy) + mean(df_01_βNS_opt.time_exact)
            list_H_hat[iβ,iN,iS]     = mean(df_PBD_βNS.H)
            list_n_submodular_cuts[iβ,iN,iS] = mean(df_PBD_βNS.n_submodular_cuts)
            list_P[iβ,iN,iS]         = mean(df_PBD_βNS.P)
            list_P1[iβ,iN,iS]        = mean(df_PBD_βNS.P1)
            list_δ[iβ,iN,iS]         = mean(df_PBD_βNS.δ)
            tlim_reached_PBD[iβ,iN,iS] = size(filter(row -> row.time_exact ≥ row.limit_time, df_PBD_βNS))[1]
            tlim_reached_SAA[iβ,iN,iS]  = size(filter(row -> row.time_exact ≥ row.limit_time, df_SAA_βNS))[1]
            tlim_reached_01[iβ,iN,iS]  = size(filter(row -> row.time_exact ≥ row.limit_time, df_01_βNS))[1]
            list_params[iβ,iN,iS]    = "(βmin,βmax,N,S) = ("*string(βmin)*","*string(βmax)*","*string(N)*","*string(S)*")"
            for ib ∈ 1:nb
                b = list_b[ib]
                df_βNSb = filter(row -> (row.α==α && row.b==b), df_βNS)
                df_PBD_βNSb = filter(row -> (row.P1_auto), df_βNSb)
                list_Z_MNL_SIM_b[iβ,iN,iS,ib] = 100*mean(df_PBD_βNSb.Z_opt_stoch)/N #objective value stored in % of capture
                list_Z_SIM_b[iβ,iN,iS,ib] = 100*mean(df_PBD_βNSb.Z_opt)/N #objective value stored in % of capture

                #estimate the value of the SIM solution for the MMNL instance
                n_rep=0
                for x_opt_str ∈ df_PBD_βNSb.x_opt
                    x_opt = eval(Meta.parse(x_opt_str))
                    list_Z_MMNL_SIM_b[iβ,iN,iS,ib] += 100*objective_value_MNL(x_opt,q,V,W_super)/N_reference #objective value stored in % of capture
                    n_rep += 1
                end
                list_Z_MMNL_SIM_b[iβ,iN,iS,ib] /= n_rep
                list_RGap_SIM_b[iβ,iN,iS,ib] = 100*(list_Z_MNL_SIM_b[iβ,iN,iS,ib]-list_Z_MMNL_SIM_b[iβ,iN,iS,ib])/list_Z_MNL_SIM_b[iβ,iN,iS,ib]
            end
            list_Z_SIM[iβ,iN,iS]      = mean(list_Z_SIM_b[iβ,iN,iS,:])
            list_Z_MNL_SIM[iβ,iN,iS]  = mean(list_Z_MNL_SIM_b[iβ,iN,iS,:])
            list_Z_MMNL_SIM[iβ,iN,iS] = mean(list_Z_MMNL_SIM_b[iβ,iN,iS,:])
            list_RGap_SIM[iβ,iN,iS]   = mean(list_RGap_SIM_b[iβ,iN,iS,:])
        end
        list_H_hat_avg[iβ] = mean(list_H_hat[iβ,:,:])
    end

    for iN ∈ 1:nN_GGX
        N = list_N_GGX[iN]
        df_GGX_βN = filter(row -> row.N==N, df_GGX_β)
        list_params_GGX[iβ,iN]    = "(βmin,βmax,N) = ("*string(βmin)*","*string(βmax)*","*string(N)*")"
        time_GGX[iβ,iN] = mean(df_GGX_βN.time)

        for ib ∈ 1:nb
            b = list_b[ib]
            df_GGX_βNb = filter(row -> (row.r==b), df_GGX_βN)
            list_Z_MNL_GGX_b[iβ,iN,ib] = 100*mean(df_GGX_βNb[!,"Z_N(X_N_GGX)"])/N #objective value stored in % of capture

            #estimate the value of the GGX solution for the MMNL instance
            n_rep=0
            for x_opt_str ∈ df_GGX_βNb.solution
                x_opt = eval(Meta.parse(x_opt_str))
                list_Z_MMNL_GGX_b[iβ,iN,ib] += 100*objective_value_MNL(x_opt,q,V,W_super)/N_reference #objective value stored in % of capture
                n_rep += 1
            end
            list_Z_MMNL_GGX_b[iβ,iN,ib] /= n_rep
            list_RGap_GGX_b[iβ,iN,ib] = 100*(list_Z_MNL_GGX_b[iβ,iN,ib]-list_Z_MMNL_GGX_b[iβ,iN,ib])/list_Z_MNL_GGX_b[iβ,iN,ib]
        end
        list_Z_MNL_GGX[iβ,iN]  = mean(list_Z_MNL_GGX_b[iβ,iN,:])
        list_Z_MMNL_GGX[iβ,iN] = mean(list_Z_MMNL_GGX_b[iβ,iN,:])
        list_RGap_GGX[iβ,iN]   = mean(list_RGap_GGX_b[iβ,iN,:])
    end

    for iN ∈ 1:nN_MOA
        N = list_N_MOA[iN]
        df_MOA_βN = filter(row -> row.N==N, df_MOA_β)
        list_params_MOA[iβ,iN] = "(βmin,βmax,N) = ("*string(βmin)*","*string(βmax)*","*string(N)*")"
        time_MOA[iβ,iN] = mean(df_MOA_βN.time)

        for ib ∈ 1:nb
            b = list_b[ib]
            df_MOA_βNb = filter(row -> (row.r==b), df_MOA_βN)
            list_Z_MNL_MOA_b[iβ,iN,ib] = 100*mean(df_MOA_βNb[!,"Z_N(X_N)"])/N #objective value stored in % of capture

            #estimate the value of the MOA solution for the MMNL instance
            n_rep=0
            for x_opt_str ∈ df_MOA_βNb.solution
                x_opt = eval(Meta.parse(x_opt_str))
                list_Z_MMNL_MOA_b[iβ,iN,ib] += 100*objective_value_MNL(x_opt,q,V,W_super)/N_reference #objective value stored in % of capture
                n_rep += 1
            end
            list_Z_MMNL_MOA_b[iβ,iN,ib] /= n_rep
            list_RGap_MOA_b[iβ,iN,ib] = 100*(list_Z_MNL_MOA_b[iβ,iN,ib]-list_Z_MMNL_MOA_b[iβ,iN,ib])/list_Z_MNL_MOA_b[iβ,iN,ib]
        end
        list_Z_MNL_MOA[iβ,iN]  = mean(list_Z_MNL_MOA_b[iβ,iN,:])
        list_Z_MMNL_MOA[iβ,iN] = mean(list_Z_MMNL_MOA_b[iβ,iN,:])
        list_RGap_MOA[iβ,iN]   = mean(list_RGap_MOA_b[iβ,iN,:])
    end
end

# Print the table content to standard output
println("____________________ Content Tables generative MMNL MIX ____________________")
println("list_params       = ",list_params)
println("list_βmin         = ",list_βmin)
println("list_βmax         = ",list_βmax)
println("list_H_hat        = ",round.(list_H_hat, digits=2))
println("list_H_hat_avg    = ",round.(list_H_hat_avg, digits=2))

println("\n................... SIM (PBD/SAA/0-1) ...................")
println("list_Z_SIM       = ",round.(list_Z_SIM, digits=2))
println("list_Z_MNL_SIM   = ",round.(list_Z_MNL_SIM, digits=2))
println("list_Z_MMNL_SIM  = ",round.(list_Z_MMNL_SIM, digits=2))
println("list_RGap_SIM    = ",round.(list_RGap_SIM, digits=2))

println("\n.......................... PBD ..........................")
println("time_exact_PBD         = ",round.(time_exact_PBD, digits=3))
println("time_total_PBD         = ",round.(time_total_PBD, digits=3))
println("tlim_reached_PBD       = ",round.(tlim_reached_PBD, digits=2))
println("list_P                 = ",Int.(round.(list_P)))
println("list_P1                = ",Int.(round.(list_P1)))
println("list_δ                 = ",round.(list_δ, digits=2))
println("list_n_submodular_cuts = ",round.(list_n_submodular_cuts, digits=1))
println("list_n_nodes_BnC_PBD   = ",round.(list_n_nodes_BnC_PBD, digits=1))

println("\n........................ SAA .......................")
println("time_exact_SAA       = ",round.(time_exact_SAA, digits=3))
println("time_total_SAA       = ",round.(time_total_SAA, digits=3))
println("tlim_reached_SAA     = ",round.(tlim_reached_SAA, digits=2))
println("list_n_nodes_BnC_SAA = ",round.(list_n_nodes_BnC_SAA, digits=1))

println("\n..................... Pure 0-1 .....................")
println("time_exact_01        = ",round.(time_exact_01, digits=3))
println("time_total_01        = ",round.(time_total_01, digits=3))
println("tlim_reached_01      = ",round.(tlim_reached_01, digits=2))
println("list_n_nodes_BnC_S01 = ",round.(list_n_nodes_BnC_01, digits=1))

println("\n........................ GGX .......................")
println("list_params_GGX  = ",list_params_GGX)
println("list_Z_MNL_GGX   = ",round.(list_Z_MNL_GGX, digits=2))
println("list_Z_MMNL_GGX  = ",round.(list_Z_MMNL_GGX, digits=2))
println("list_RGap_GGX    = ",round.(list_RGap_GGX, digits=2))
println("time_GGX         = ",round.(time_GGX, digits=2))

println("\n........................ MOA .......................")
println("list_params_MOA  = ",list_params_MOA)
println("list_Z_MNL_MOA   = ",round.(list_Z_MNL_MOA, digits=2))
println("list_Z_MMNL_MOA  = ",round.(list_Z_MMNL_MOA, digits=2))
println("list_RGap_MOA    = ",round.(list_RGap_MOA, digits=2))
println("time_MOA         = ",round.(time_MOA, digits=2))