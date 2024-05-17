# Author: Robin Legault <legault@mit.edu>

# Packages
using DataFrames, XLSX, Statistics, Plots, LaTeXStrings, Random

# Dependencies
include("../src/IO_functions.jl")

##########################################################################################################
#################### Computational results for conditional MNL instances, Section 6.1. ###################
##########################################################################################################

results_NYC = true
results_HM14 = true

if(results_NYC)
    # Load results
    results_file = "MNL_NYC_paper.xlsx" #
    df = DataFrame(XLSX.readtable(RESULTS_PATH*results_file, "Sheet1"; header=true, infer_eltypes=true))
    df_MOA = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*"GGX_MOA_"*results_file, "MOA"; header=true, infer_eltypes=true))
    df_GGX = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*"GGX_MOA_"*results_file, "GGX"; header=true, infer_eltypes=true))
    col_β = parse.(Float64, [match(r"beta_([\d.]+)", string).captures[1] for string in df.inst_name])
    col_α = parse.(Float64, [match(r"alpha_([\d.]+)", string).captures[1] for string in df.inst_name])
    df[!, "β"] = col_β
    df[!, "α"] = col_α
    N = df.N[1]
    list_β = reverse(sort!(unique(df.β)))
    list_S = unique(df.S)
    list_α = unique(df.α)
    list_b = unique(df.b)
    nβ = length(list_β)
    nS = length(list_S)
    nα = length(list_α)
    nb = length(list_b)

    # Collect statistics
    time_exact_PBD         = zeros(nβ, nS)
    time_exact_SAA         = zeros(nβ, nS)
    time_exact_01          = zeros(nβ, nS)
    time_total_PBD         = zeros(nβ, nS)
    time_total_SAA         = zeros(nβ, nS)
    time_total_01          = zeros(nβ, nS)
    tlim_reached_PBD       = zeros(nβ, nS)
    tlim_reached_SAA       = zeros(nβ, nS)
    tlim_reached_01        = zeros(nβ, nS)
    list_H_hat             = zeros(nβ, nS)
    list_H_hat_avg         = zeros(nβ)
    list_n_submodular_cuts = zeros(nβ, nS)
    list_n_nodes_BnC_PBD   = zeros(nβ, nS)
    list_n_nodes_BnC_SAA   = zeros(nβ, nS)
    list_n_nodes_BnC_01    = zeros(nβ, nS)
    list_P                 = zeros(nβ, nS)
    list_P1                = zeros(nβ, nS)
    list_δ                 = zeros(nβ, nS)
    list_RGap_GGX          = zeros(nβ)
    list_RGap_SIM          = zeros(nβ, nS)
    list_Z_MNL_MOA_αb      = zeros(nβ, nα, nb)
    list_Z_MNL_GGX_αb      = zeros(nβ, nα, nb)
    list_RGap_GGX_αb       = zeros(nβ, nα, nb)
    list_Z_SIM_αb          = zeros(nβ, nS, nα, nb)
    list_Z_MNL_SIM_αb      = zeros(nβ, nS, nα, nb)
    list_RGap_SIM_αb       = zeros(nβ, nS, nα, nb)
    list_time_MOA          = zeros(nβ)
    list_time_GGX          = zeros(nβ)
    list_params            = Array{String}(undef, nβ, nS)
    for iβ ∈ 1:nβ
        β = list_β[iβ]
        df_β = filter(row -> row.β==β, df)
        df_MOA_β = filter(row -> row.beta==β, df_MOA)
        df_GGX_β = filter(row -> row.beta==β, df_GGX)
        list_time_MOA[iβ] = mean(df_MOA_β[!,"time"])
        list_time_GGX[iβ] = mean(df_GGX_β[!,"time"])
        for iα ∈ 1:nα
            α = list_α[iα]
            for ib ∈ 1:nb
                b = list_b[ib]
                df_MOA_βαb = filter(row -> (row.alpha==α && row.r==b), df_MOA_β)
                df_GGX_βαb = filter(row -> (row.alpha==α && row.r==b), df_GGX_β)
                list_Z_MNL_MOA_αb[iβ,iα,ib] = df_MOA_βαb[!,"Z_N(X_N)"][1]
                list_Z_MNL_GGX_αb[iβ,iα,ib] = df_GGX_βαb[!,"Z_N(X_N_GGX3)"][1]
                list_RGap_GGX_αb[iβ,iα,ib] = 100*(list_Z_MNL_MOA_αb[iβ,iα,ib]-list_Z_MNL_GGX_αb[iβ,iα,ib])/list_Z_MNL_MOA_αb[iβ,iα,ib]
            end
        end
        list_RGap_GGX[iβ] = mean(list_RGap_GGX_αb[iβ,:,:])
        for iS ∈ 1:nS
            S = list_S[iS]
            df_βS  = filter(row -> row.S==S, df_β)
            df_PBD_βS = filter(row -> (row.P1_auto), df_βS)
            df_SAA_βS  = filter(row -> (row.time_clustering==0.0), df_βS)
            df_01_βS  = filter(row -> (row.time_clustering>0.0 && !row.P1_auto && row.P1_weight_min==1), df_βS)
            df_PBD_βS_opt = filter(row -> row.is_optimal, df_PBD_βS)
            df_SAA_βS_opt  = filter(row -> row.is_optimal, df_SAA_βS)
            df_01_βS_opt  = filter(row -> row.is_optimal, df_01_βS)
            list_n_nodes_BnC_PBD[iβ,iS]   = mean(df_PBD_βS_opt.n_nodes_BnC)
            list_n_nodes_BnC_SAA[iβ,iS]   = mean(df_SAA_βS_opt.n_nodes_BnC)
            list_n_nodes_BnC_01[iβ,iS]    = mean(df_01_βS_opt.n_nodes_BnC)
            time_exact_PBD[iβ,iS] = mean(df_PBD_βS_opt.time_exact)
            time_exact_SAA[iβ,iS]  = mean(df_SAA_βS_opt.time_exact)
            time_exact_01[iβ,iS]  = mean(df_01_βS_opt.time_exact)
            time_total_PBD[iβ,iS] = mean(df_PBD_βS_opt.time_clustering) + mean(df_PBD_βS_opt.time_greedy) + mean(df_PBD_βS_opt.time_exact)
            time_total_SAA[iβ,iS]  = mean(df_SAA_βS_opt.time_clustering)  + mean(df_SAA_βS_opt.time_greedy)  + mean(df_SAA_βS_opt.time_exact)
            time_total_01[iβ,iS]  = mean(df_01_βS_opt.time_clustering)  + mean(df_01_βS_opt.time_greedy)  + mean(df_01_βS_opt.time_exact)
            list_H_hat[iβ,iS]     = mean(df_PBD_βS.H)
            list_n_submodular_cuts[iβ,iS] = mean(df_PBD_βS.n_submodular_cuts)
            list_P[iβ,iS]         = mean(df_PBD_βS.P)
            list_P1[iβ,iS]        = mean(df_PBD_βS.P1)
            list_δ[iβ,iS]         = mean(df_PBD_βS.δ)
            tlim_reached_PBD[iβ,iS] = size(filter(row -> row.time_exact ≥ row.limit_time, df_PBD_βS))[1]
            tlim_reached_SAA[iβ,iS]  = size(filter(row -> row.time_exact ≥ row.limit_time, df_SAA_βS))[1]
            tlim_reached_01[iβ,iS]  = size(filter(row -> row.time_exact ≥ row.limit_time, df_01_βS))[1]
            list_params[iβ,iS]    = "(β,S) = ("*string(β)*","*string(S)*")"
            for iα ∈ 1:nα
                α = list_α[iα]
                for ib ∈ 1:nb
                    b = list_b[ib]
                    df_βSαb = filter(row -> (row.α==α && row.b==b), df_βS)
                    df_PBD_βSαb = filter(row -> (row.P1_auto), df_βSαb)
                    list_Z_SIM_αb[iβ,iS,iα,ib]     = mean(df_PBD_βSαb.Z_opt)
                    list_Z_MNL_SIM_αb[iβ,iS,iα,ib] = mean(df_PBD_βSαb.Z_opt_stoch)
                    list_RGap_SIM_αb[iβ,iS,iα,ib]  = 100*(list_Z_MNL_MOA_αb[iβ,iα,ib]-list_Z_MNL_SIM_αb[iβ,iS,iα,ib])/list_Z_MNL_MOA_αb[iβ,iα,ib]
                end
            end
            list_RGap_SIM[iβ,iS] = mean(list_RGap_SIM_αb[iβ,iS,:,:])
        end
        list_H_hat_avg[iβ] = mean(list_H_hat[iβ,:])
    end


    # Print the table content to standard output
    println("____________________ Content Tables conditional MNL NYC ____________________")
    println("list_params    = ",list_params)
    println("list_β         = ",list_β)
    println("list_H_hat     = ",round.(list_H_hat, digits=2))
    println("list_H_hat_avg = ",round.(list_H_hat_avg, digits=2))

    println("\n................... SIM (PBD/SAA/0-1) ...................")
    println("list_RGap_SIM  = ",round.(list_RGap_SIM, digits=2))

    println("\n.......................... PBD ..........................")
    println("time_exact_PBD = ",round.(time_exact_PBD, digits=2))
    println("time_total_PBD = ",round.(time_total_PBD, digits=2))
    println("tlim_reached_PBD = ",round.(tlim_reached_PBD, digits=2))
    println("list_P         = ",Int.(round.(list_P)))
    println("list_P1        = ",Int.(round.(list_P1)))
    println("100*P/R        = ",round.(100*list_P./(N*list_S)', digits=1))
    println("100*P1/P       = ",round.(100*list_P1./list_P, digits=1))
    println("list_δ(%)      = ",round.(100*list_δ, digits=1))
    println("list_n_submodular_cuts = ",round.(list_n_submodular_cuts, digits=1))
    println("list_n_nodes_BnC_PBD = ",round.(list_n_nodes_BnC_PBD, digits=1))

    println("\n........................ SAA .......................")
    println("time_exact_SAA       = ",round.(time_exact_SAA, digits=2))
    println("time_total_SAA       = ",round.(time_total_SAA, digits=2))
    println("tlim_reached_SAA     = ",round.(tlim_reached_SAA, digits=2))
    println("list_n_nodes_BnC_SAA = ",round.(list_n_nodes_BnC_SAA, digits=1))

    println("\n..................... Pure 0-1 .....................")
    println("time_exact_01        = ",round.(time_exact_01, digits=2))
    println("time_total_01        = ",round.(time_total_01, digits=2))
    println("tlim_reached_01      = ",round.(tlim_reached_01, digits=2))
    println("list_n_nodes_BnC_S01 = ",round.(list_n_nodes_BnC_01, digits=1))

    println("\n........................ MOA .......................")
    println("list_time_MOA  = ",round.(list_time_MOA, digits=2))

    println("\n........................ GGX .......................")
    println("list_time_GGX  = ",round.(list_time_GGX, digits=2))
    println("list_RGap_GGX  = ",round.(list_RGap_GGX , digits=2))
end

if(results_HM14)
    # Load results
    results_file = "MNL_HM14_paper.xlsx" 
    df = DataFrame(XLSX.readtable(RESULTS_PATH*results_file, "Sheet1"; header=true, infer_eltypes=true))
    df_MOA = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*"GGX_MOA_"*results_file, "MOA"; header=true, infer_eltypes=true))
    df_GGX = DataFrame(XLSX.readtable(RESULTS_PATH_MOA_GGX*"GGX_MOA_"*results_file, "GGX"; header=true, infer_eltypes=true))
    col_β = parse.(Float64, [match(r"beta_([\d.]+)", string).captures[1] for string in df.inst_name])
    col_α = parse.(Float64, [match(r"alpha_([\d.]+)", string).captures[1] for string in df.inst_name])
    df[!, "β"] = col_β
    df[!, "α"] = col_α
    N = df.N[1]
    list_D = unique(df.D)
    list_β = reverse(sort!(unique(df.β)))
    list_S = unique(df.S)
    list_α = unique(df.α)
    list_b = unique(df.b)
    nD = length(list_D)
    nβ = length(list_β)
    nS = length(list_S)
    nα = length(list_α)
    nb = length(list_b)

    # Collect statistics
    time_exact_PBD         = zeros(nD, nβ, nS)
    time_exact_SAA         = zeros(nD, nβ, nS)
    time_exact_01          = zeros(nD, nβ, nS)
    time_total_PBD         = zeros(nD, nβ, nS)
    time_total_SAA         = zeros(nD, nβ, nS)
    time_total_01          = zeros(nD, nβ, nS)
    tlim_reached_PBD       = zeros(nD, nβ, nS)
    tlim_reached_SAA       = zeros(nD, nβ, nS)
    tlim_reached_01        = zeros(nD, nβ, nS)
    list_H_hat             = zeros(nD, nβ, nS)
    list_H_hat_avg         = zeros(nD, nβ)
    list_n_submodular_cuts = zeros(nD, nβ, nS)
    list_n_nodes_BnC_PBD   = zeros(nD, nβ, nS)
    list_n_nodes_BnC_SAA   = zeros(nD, nβ, nS)
    list_n_nodes_BnC_01    = zeros(nD, nβ, nS)
    list_P                 = zeros(nD, nβ, nS)
    list_P1                = zeros(nD, nβ, nS)
    list_δ                 = zeros(nD, nβ, nS)
    list_RGap_GGX          = zeros(nD, nβ)
    list_RGap_SIM          = zeros(nD, nβ, nS)
    list_Z_MNL_best_αb     = zeros(nD, nβ, nα, nb)
    list_Z_MNL_MOA_αb      = zeros(nD, nβ, nα, nb)
    list_Z_MNL_GGX_αb      = zeros(nD, nβ, nα, nb)
    list_RGap_GGX_αb       = zeros(nD, nβ, nα, nb)
    list_Z_SIM_αb          = zeros(nD, nβ, nS, nα, nb)
    list_Z_MNL_SIM_αb      = zeros(nD, nβ, nS, nα, nb)
    list_RGap_SIM_αb       = zeros(nD, nβ, nS, nα, nb)
    tlim_reached_01_αb     = zeros(Bool, nD, nβ, nS, nα, nb)
    list_time_MOA          = zeros(nD, nβ)
    list_time_GGX          = zeros(nD, nβ)
    list_params            = Array{String}(undef, nD, nβ, nS)
    for iD ∈ 1:nD
        D = list_D[iD]
        df_D = filter(row -> row.D==D, df)
        df_MOA_D = filter(row -> row.D==D, df_MOA)
        df_GGX_D = filter(row -> row.D==D, df_GGX)
        for iβ ∈ 1:nβ
            β = list_β[iβ]
            df_Dβ = filter(row -> row.β==β, df_D)
            df_MOA_Dβ = filter(row -> row.beta==β, df_MOA_D)
            df_GGX_Dβ = filter(row -> row.beta==β, df_GGX_D)
            list_time_MOA[iD,iβ] = mean(df_MOA_Dβ[!,"time"])
            list_time_GGX[iD,iβ] = mean(df_GGX_Dβ[!,"time"])
            for iα ∈ 1:nα
                α = list_α[iα]
                for ib ∈ 1:nb
                    b = list_b[ib]
                    df_MOA_Dβαb = filter(row -> (row.alpha==α && row.r==b), df_MOA_Dβ)
                    df_GGX_Dβαb = filter(row -> (row.alpha==α && row.r==b), df_GGX_Dβ)
                    list_Z_MNL_MOA_αb[iD,iβ,iα,ib] = df_MOA_Dβαb[!,"Z_N(X_N)"][1]
                    list_Z_MNL_GGX_αb[iD,iβ,iα,ib] = df_GGX_Dβαb[!,"Z_N(X_N_GGX3)"][1]
                end
            end

            for iS ∈ 1:nS
                S = list_S[iS]
                df_DβS  = filter(row -> row.S==S, df_Dβ)
                df_PBD_DβS = filter(row -> (row.P1_auto), df_DβS)
                df_SAA_DβS  = filter(row -> (row.time_clustering==0.0), df_DβS)
                df_01_DβS  = filter(row -> (row.time_clustering>0.0 && !row.P1_auto && row.P1_weight_min==1), df_DβS)
                df_PBD_DβS_opt = filter(row -> row.is_optimal, df_PBD_DβS)
                df_SAA_DβS_opt  = filter(row -> row.is_optimal, df_SAA_DβS)
                df_01_DβS_opt  = filter(row -> row.is_optimal, df_01_DβS)
                list_n_nodes_BnC_PBD[iD,iβ,iS] = mean(df_PBD_DβS_opt.n_nodes_BnC)
                list_n_nodes_BnC_SAA[iD,iβ,iS] = mean(df_SAA_DβS_opt.n_nodes_BnC)
                list_n_nodes_BnC_01[iD,iβ,iS]  = mean(df_01_DβS_opt.n_nodes_BnC)
                time_exact_PBD[iD,iβ,iS] = mean(df_PBD_DβS_opt.time_exact)
                time_exact_SAA[iD,iβ,iS]  = mean(df_SAA_DβS_opt.time_exact)
                time_exact_01[iD,iβ,iS]  = mean(df_01_DβS_opt.time_exact)
                time_total_PBD[iD,iβ,iS] = mean(df_PBD_DβS_opt.time_clustering) + mean(df_PBD_DβS_opt.time_greedy) + mean(df_PBD_DβS_opt.time_exact)
                time_total_SAA[iD,iβ,iS]  = mean(df_SAA_DβS_opt.time_clustering)  + mean(df_SAA_DβS_opt.time_greedy)  + mean(df_SAA_DβS_opt.time_exact)
                time_total_01[iD,iβ,iS]  = mean(df_01_DβS_opt.time_clustering)  + mean(df_01_DβS_opt.time_greedy)  + mean(df_01_DβS_opt.time_exact)
                list_H_hat[iD,iβ,iS]     = mean(df_PBD_DβS.H)
                list_n_submodular_cuts[iD,iβ,iS] = mean(df_PBD_DβS.n_submodular_cuts)
                list_P[iD,iβ,iS]         = mean(df_PBD_DβS.P)
                list_P1[iD,iβ,iS]        = mean(df_PBD_DβS.P1)
                list_δ[iD,iβ,iS]         = mean(df_PBD_DβS.δ)
                tlim_reached_PBD[iD,iβ,iS] = size(filter(row -> row.time_exact ≥ row.limit_time, df_PBD_DβS))[1]
                tlim_reached_SAA[iD,iβ,iS]  = size(filter(row -> row.time_exact ≥ row.limit_time, df_SAA_DβS))[1]
                tlim_reached_01[iD,iβ,iS]  = size(filter(row -> row.time_exact ≥ row.limit_time, df_01_DβS))[1]
                list_params[iD,iβ,iS]    = "(D,β,S) = ("*string(D)*","*string(β)*","*string(S)*")"
                for iα ∈ 1:nα
                    α = list_α[iα]
                    for ib ∈ 1:nb
                        b = list_b[ib]
                        df_DβSαb = filter(row -> (row.α==α && row.b==b), df_DβS)
                        df_01_DβSαb = filter(row -> (row.time_clustering>0.0 && !row.P1_auto && row.P1_weight_min==1), df_DβSαb)
                        list_Z_SIM_αb[iD,iβ,iS,iα,ib]     = mean(df_01_DβSαb.Z_opt)
                        list_Z_MNL_SIM_αb[iD,iβ,iS,iα,ib] = mean(df_01_DβSαb.Z_opt_stoch)
                        tlim_reached_01_αb[iD,iβ,iS,iα,ib] = !df_01_DβSαb.is_optimal[1]
                    end
                end
            end
            list_H_hat_avg[iD,iβ] = mean(list_H_hat[iD,iβ,:])

            for iα ∈ 1:nα
                for ib ∈ 1:nb
                    #take the best known value for each instance (list_Z_MNL_MOA_αb[iD,iβ,iα,ib] in theory, but MOA has numerical stability issues that can lead to suboptimal solutions)
                    list_Z_MNL_best_αb[iD,iβ,iα,ib] = max(list_Z_MNL_MOA_αb[iD,iβ,iα,ib], list_Z_MNL_GGX_αb[iD,iβ,iα,ib], maximum(list_Z_MNL_SIM_αb[iD,iβ,:,iα,ib]))
                    list_RGap_GGX_αb[iD,iβ,iα,ib] = 100*(list_Z_MNL_best_αb[iD,iβ,iα,ib]-list_Z_MNL_GGX_αb[iD,iβ,iα,ib])/list_Z_MNL_best_αb[iD,iβ,iα,ib]
                    for iS ∈ 1:nS
                        list_RGap_SIM_αb[iD,iβ,iS,iα,ib] = 100*(list_Z_MNL_best_αb[iD,iβ,iα,ib]-list_Z_MNL_SIM_αb[iD,iβ,iS,iα,ib])/list_Z_MNL_best_αb[iD,iβ,iα,ib]
                    end
                end
            end
            list_RGap_GGX[iD,iβ] = mean(list_RGap_GGX_αb[iD,iβ,:,:])

            #RGAP for simulation-based methods only computed based on instances that could be solved within the time limit
            for iS ∈ 1:nS
                sum_RGap = 0.0
                count_RGap = 0
                for iα ∈ 1:nα
                    for ib ∈ 1:nb
                        if(!tlim_reached_01_αb[iD,iβ,iS,iα,ib])
                            sum_RGap += list_RGap_SIM_αb[iD,iβ,iS,iα,ib]
                            count_RGap += 1
                        end
                    end
                end
                list_RGap_SIM[iD,iβ,iS] = sum_RGap/count_RGap
            end
        end
    end

    # Print the table content to standard output
    println("____________________ Content Tables conditional MNL HM14 ____________________")
    for iD ∈ 1:nD
        D = list_D[iD]
        println("\n========================== D = ",D," ==========================")
        println("list_params    = ",list_params[iD,:,:])
        println("list_β         = ",list_β)
        println("list_S         = ",list_S)
        println("list_H_hat     = ",round.(list_H_hat[iD,:,:], digits=2))
        println("list_H_hat_avg = ",round.(list_H_hat_avg[iD,:], digits=2))

        println("\n................... SIM (PBD/SAA/0-1) ...................")
        println("list_RGap_SIM  = ",round.(list_RGap_SIM[iD,:,:], digits=2))
    
        println("\n.......................... PBD ..........................")
        println("time_exact_PBD         = ",round.(time_exact_PBD[iD,:,:], digits=3))
        println("time_total_PBD         = ",round.(time_total_PBD[iD,:,:], digits=3))
        println("tlim_reached_PBD       = ",round.(tlim_reached_PBD[iD,:,:]))
        println("list_P                 = ",Int.(round.(list_P[iD,:,:])))
        println("list_P1                = ",Int.(round.(list_P1[iD,:,:])))
        println("100*P/R                = ",round.(100*list_P[iD,:,:]./(N*list_S)', digits=1))
        println("100*P1/P               = ",round.(100*list_P1[iD,:,:]./list_P[iD,:,:], digits=1))
        println("list_δ(%)              = ",round.(100*list_δ[iD,:,:], digits=1))
        println("list_n_submodular_cuts = ",round.(list_n_submodular_cuts[iD,:,:], digits=1))
        println("list_n_nodes_BnC_PBD   = ",round.(list_n_nodes_BnC_PBD[iD,:,:], digits=1))


        println("\n........................ SAA .......................")
        println("time_exact_SAA       = ",round.(time_exact_SAA[iD,:,:], digits=3))
        println("time_total_SAA       = ",round.(time_total_SAA[iD,:,:], digits=3))
        println("tlim_reached_SAA     = ",round.(tlim_reached_SAA[iD,:,:]))
        println("list_n_nodes_BnC_SAA = ",round.(list_n_nodes_BnC_SAA[iD,:,:], digits=1))

        println("\n..................... Pure 0-1 .....................")
        println("time_exact_01       = ",round.(time_exact_01[iD,:,:], digits=3))
        println("time_total_01       = ",round.(time_total_01[iD,:,:], digits=3))
        println("tlim_reached_01     = ",round.(tlim_reached_01[iD,:,:]))
        println("list_n_nodes_BnC_01 = ",round.(list_n_nodes_BnC_01[iD,:,:], digits=1))

        println("\n........................ MOA .......................")
        println("list_time_MOA  = ",round.(list_time_MOA[iD,:], digits=2))

        println("\n........................ GGX .......................")
        println("list_time_GGX  = ",round.(list_time_GGX[iD,:], digits=2))
        println("list_RGap_GGX  = ",round.(list_RGap_GGX[iD,:], digits=2))
    end
end