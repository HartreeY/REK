using StatsBase, Distributions, Random, Plots, SpecialFunctions, Serialization, Dates

# Constants
# ------------------------------------------------

# Plot configuration
default(margin=6Plots.mm)

# Space parameters
DEF_X_MAX = 100
DEF_Y_MAX = 10
DEF_Z_MAX = 10

# Population parameters
DEF_N_DEMES_STARTFILL = 5
DEF_K_CAPACITY = 20
DEF_R_PROLIF_RATE = 1.8
DEF_r_LOG_PROLIF_RATE = log(DEF_R_PROLIF_RATE)

# Gene parameters (regular)
DEF_N_LOCI = 1000
DEF_N_SEL_LOCI = 500
DEF_MUT_RATE = 0.7567 # genome-wide
DEF_MIGR_RATE = 0.1
DEF_S_SEL_COEF = 0.002
DEF_H_DOMIN_COEF = 0
DEF_PROP_OF_DEL_MUTS = 0.9
MIGR_DIRS_2D_8 = [[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]
MIGR_DIRS_2D_6 = [[-1,0],[0,-1],[-1,1],[0,1],[1,0],[1,1]]
MIGR_DIRS_2D_4 = [[-1,0],[0,-1],[0,1],[1,0]]
MIGR_DIRS_3D_8 = [[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]
MIGR_DIRS_3D_6 = [[-1,0],[0,-1],[-1,1],[0,1],[1,0],[1,1]]
MIGR_DIRS_3D_4 = [[-1,0],[0,-1],[0,1],[1,0]]

# Gene parameters (infinite-sites)
DEF_N_SEGR_REGIONS = 20
DEF_PROP_OF_SEL_LOCI = 1.0

# Expansion parameters
DEF_X_MAX_BURNIN = 5
DEF_N_GENS_BURNIN = 10
DEF_X_MAX_EXP = DEF_X_MAX
DEF_Y_MAX_EXP = DEF_Y_MAX
DEF_N_GENS_EXP = 20
DEF_MIGR_MODE_2D = "4"
DEF_MIGR_MODE_3D = "6"
DEF_DATA_TO_GENERATE = "FP"

# Structures
# ------------------------------------------------

"""
Calculates average fitness and average hetero- and homozygosities in a deme in the finite-site model.
"""
function calc_muts_and_meanf_in_deme(ms1,ms2,s_sel_coef,h_domin_coef,n_loci,sel_loci=[])
    len = length(ms1)
    muts_AAsel_total = 0
    muts_Aasel_total = 0
    muts_aasel_total = 0
    muts_AAneu_total = 0
    muts_Aaneu_total = 0
    muts_aaneu_total = 0
    fits = []
    
    for i in 1:len
        muts_AA_sel = 0
        muts_Aa_sel = 0
        muts_AA_neu = 0
        muts_Aa_neu = 0
        new_fitness = 1.0

        for j in 1:n_loci
            if ms1[i][j]==true && ms2[i][j]==true
                if j in sel_loci
                    muts_AA_sel += 1
                    new_fitness *= 1 - s_sel_coef
                else
                    muts_AA_neu += 1
                end

            elseif ms1[i][j]==true || ms2[i][j]==true
                if j in sel_loci
                    muts_Aa_sel += 1
                    new_fitness *= 1 - h_domin_coef * s_sel_coef
                else
                    muts_Aa_neu += 1
                end
            end
        end

        push!(fits,new_fitness)
        muts_AAsel_total += muts_AA_sel
        muts_Aasel_total += muts_Aa_sel
        muts_AAneu_total += muts_AA_neu
        muts_Aaneu_total += muts_Aa_neu
        muts_aasel_total += length(sel_loci) - muts_AA_sel - muts_Aa_sel
        muts_aaneu_total += n_loci - length(sel_loci) - muts_AA_neu - muts_Aa_neu
    end

    muts_AAsel_total /= len
    muts_Aasel_total /= len
    muts_aasel_total /= len
    muts_AAneu_total /= len
    muts_Aaneu_total /= len
    muts_aaneu_total /= len
    return muts_AAsel_total,muts_Aasel_total,muts_aasel_total,muts_AAneu_total,muts_Aaneu_total,muts_aaneu_total,fits
end

@inbounds function mutate(ms1,ms2,mut_rate,n_loci)
    get_mutation_random = rand(Poisson(mut_rate))
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:n_loci)

        if rand(1:2)==1
            ms1[pos_alter] = true
        else
            ms2[pos_alter] = true
        end
    end
end

@inbounds function crossover(ms1,ms2,n_loci)
    for j in 1:n_loci
        lr = rand(1:2)
        ms1[j] = lr==1 ? ms1[j] : ms2[j]
    end
end

@inbounds function mate(ind1,ind2,n_loci)
    new_loci = vcat(ind1[1:n_loci],ind2[1:n_loci])
    return new_loci
end

@inbounds function calc_offspring(pnt_wld,pnt_wld_stats)
    next_gen_posits = []
    next_gen_pops = fill(NaN,pnt_wld_stats["max"]...)
    for k in Iterators.product([1:n for n in pnt_wld_stats["max"]]...)
        if isassigned(pnt_wld,k...) && length(pnt_wld[k...])>0
            n_ppl_at_deme = length(pnt_wld[k...])
            expected_offspring = n_ppl_at_deme * (pnt_wld_stats["r_prolif_rate"]/(1 + (n_ppl_at_deme*(pnt_wld_stats["r_prolif_rate"]-1))/pnt_wld_stats["k_capacity"]))
            next_gen_pops[k...] =  rand(Poisson(expected_offspring))
            if next_gen_pops[k...]>0
                push!(next_gen_posits,[k...])
            end
        end
    end
    return next_gen_posits, next_gen_pops
end