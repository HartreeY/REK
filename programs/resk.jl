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
MIGR_PROBS = [
    Dict(["ort"=>(1,0)]), #1D
    Dict(["ort"=>(1,0),"all"=>(1/2,1/2),"buffon1"=>(2/pi,1/pi),"buffon2"=>(4/3/pi,1/3/pi),"buffon3"=>(0.4244132,0.21221),"diag1/2"=>(2/3,1/3)]), # 2D. To add "hex"!
    Dict(["ort"=>(1,0),"all"=>(1/2,1/2),"buffon1"=>(2/pi,1/pi),"buffon2"=>(4/3/pi,1/3/pi),"buffon3"=>(0.4244132,0.21221),"diag1/2"=>(2/3,1/3)]) # 3D. To add "hex"! To confirm Buffon for 3d!
]
MIGR_DIRS_ORT = [
    [[1],[-1]], # 1D
    [[-1,0],[0,-1],[0,1],[1,0]], #2D
    [[-1,0,0],[1,0,0],[0,1,0],[0,-1,0],[0,0,-1],[0,0,1]] #3D
]
MIGR_DIRS_DIAG = [
    [], #1D
    [[-1,-1],[-1,1],[1,-1],[1,1]], #2D
    [[-1,-1,-1],[-1,-1,0],[-1,-1,1],[-1,0,-1],[-1,0,1],[-1,1,-1],[-1,1,0],[-1,1,1],[0,-1,-1],[0,-1,0],[0,-1,1],[0,0,0],[0,1,-1],[0,1,1],[1,-1,-1],[1,-1,0],[1,-1,1],[1,0,-1],[1,0,1],[1,1,-1],[1,1,0],[1,1,1]] #3D
]
MIGR_DIRS_HEX = [
    1=> [], #1D
    2=> [[-1,0],[0,-1],[-1,1],[0,1],[1,0],[1,1]], #2D
    3=> [] #3D
] # To do!

# Gene parameters (infinite-sites)
DEF_N_SEGR_REGIONS = 20
DEF_PROP_OF_SEL_LOCI = 1.0

# Expansion parameters
DEF_X_MAX_BURNIN = 5
DEF_R_MAX_BURNIN = 3
DEF_N_GENS_BURNIN = 10
DEF_X_MAX_EXP = DEF_X_MAX
DEF_Y_MAX_EXP = DEF_Y_MAX
DEF_R_MAX_EXP = 20
DEF_N_GENS_EXP = 40
DEF_MIGR_MODE = "ort"
DEF_DATA_TO_GENERATE = "FP"

# Other
DEF_GRAPH_SIZE = (640,640)

# Simulation functions (regular)
# ------------------------------------------------


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

@inbounds function calc_migr_dist(deme,pnt_wld_stats,migr_mode,bottleneck,max_migr=pnt_wld_stats["max"],refl_walls=false,r_max_migr=0)

    wlddim = pnt_wld_stats["wlddim"]
    max = pnt_wld_stats["max"]
    move = zeros(Int32,wlddim)
    migr_res = rand()
    if !(migr_mode in keys(MIGR_PROBS[wlddim]))
        migr_mode = "ort"
    end
    p_lat,p_diag = MIGR_PROBS[wlddim][migr_mode]

    if rand() < pnt_wld_stats["migr_rate"] && migr_res < p_lat + p_diag
        if migr_res < p_lat
            dir = copy(sample(MIGR_DIRS_ORT[wlddim]))
        elseif migr_res < p_lat + p_diag
            dir = copy(sample(MIGR_DIRS_DIAG[wlddim]))
        end

        # Raw migration results
        move = copy(dir)

        # Nullify migration on certain conditions
        # Inside certain radius check:
        if r_max_migr>0
            r_arr = [deme[i]-(max[i]-1)/2+move[i] for i in 1:wlddim]
            r2 = sum(r_arr.^2)
            if r2 > r_max_migr*r_max_migr
                #factor = r_max_migr*r_max_migr/r2
                move = [0,0] # Do [trunc(Int16, factor * move[i]) for i in 1:wlddim] in the future (with multiple-deme jumps)
            end
        end

        # Inside certain square check:
        if isa(max_migr,Tuple)

            for i in 1:wlddim
                try_move = deme[i]+move[i]
                if try_move > max_migr[i] || try_move < 1
                    move[i] = refl_walls ? -move[i] : 0
                end
            end
        end

        # Bottleneck barrier check:
        if isa(bottleneck,Tuple) && isa(bottleneck[2],Int) && bottleneck[2]>0
            if bottleneck[1]=="midhole at x="
                common_cond = deme[1]+move[1]==bottleneck[2] && deme[2]+move[2]!=ceil(max[2]/2)
                if wlddim==2 && common_cond
                    move .= 0
                elseif wlddim==3 && common_cond && deme[3]+move[3]!=ceil(max[3]/2)
                    move .= 0
                end
            elseif bottleneck[1]=="midhole at y="
                common_cond = deme[2]+move[2]==bottleneck[2] && deme[1]+move[1]!=ceil(max[1]/2)
                if wlddim==2 && common_cond
                    move .= 0
                elseif wlddim==3 && common_cond && deme[3]+move[3]!=ceil(max[3]/2)
                    move .= 0
                end
            elseif bottleneck[1]=="midhole at z=" && wlddim==3 && deme[3]+move[3]==bottleneck[2] && deme[1]+move[1]!=ceil(max[1]/2) && deme[2]+move[2]!=ceil(max[2]/2)
                move .= 0
            end
        end
    end
    return move
end

"""
Calculates average fitness and average hetero- and homozygosities in a deme in the finite-site model.
"""
function calc_muts_and_fitn_in_deme(ms1,ms2,s_sel_coef,h_domin_coef,n_loci,sel_loci=[])
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

@inbounds function build_next_gen(pnt_wld_ms1,pnt_wld_ms2,pnt_wld_stats,pnt_fitn_wld=NaN,pnt_pops_wld=NaN,pnt_muts_AAsel_wld=NaN,pnt_muts_Aasel_wld=NaN,
    pnt_muts_aasel_wld=NaN,pnt_muts_AAneu_wld=NaN,pnt_muts_Aaneu_wld=NaN,pnt_muts_aaneu_wld=NaN;
    max_migr=NaN,migr_mode=DEF_MIGR_MODE,bottleneck=NaN,refl_walls=false,r_max_migr=0)

    wlddim = pnt_wld_stats["wlddim"]

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(pnt_wld_ms1,pnt_wld_stats)
    
    # Define the world (as an array [=demes] of arrays [=individs] of two Bool arrays [=monosomes]) and the data arrays in the next generation
    wld_ms1_next = Array{Array{Array{Bool}}, wlddim}(undef,pnt_wld_stats["max"]...)
    wld_ms2_next = Array{Array{Array{Bool}}, wlddim}(undef,pnt_wld_stats["max"]...)
    mean_fitn_next = NaN
    pops_next = NaN
    muts_AAsel_next = NaN
    muts_Aasel_next = NaN
    muts_aasel_next = NaN
    muts_AAneu_next = NaN
    muts_Aaneu_next = NaN
    muts_aaneu_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    fitn_out = false
    pops_out = false
    sel_out = false
    neu_out = false

    if pnt_fitn_wld isa Array{Float32, wlddim+1}
        fitn_out = true
        mean_fitn_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(mean_fitn_next,NaN)
    end
    if pnt_pops_wld isa Array{Float32, wlddim+1}
        pops_out = true
        pops_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(pops_next,NaN)
    end
    if (pnt_muts_AAsel_wld isa Array{Float32, wlddim+1}) && (pnt_muts_Aasel_wld isa Array{Float32, wlddim+1}) && (pnt_muts_aasel_wld isa Array{Float32, wlddim+1})
        sel_out = true
        muts_AAsel_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        muts_Aasel_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        muts_aasel_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(muts_AAsel_next,NaN)
        fill!(muts_Aasel_next,NaN)
        fill!(muts_aasel_next,NaN)
    end
    if (pnt_muts_AAneu_wld isa Array{Float32, wlddim+1}) && (pnt_muts_Aaneu_wld isa Array{Float32, wlddim+1}) && (pnt_muts_aaneu_wld isa Array{Float32, wlddim+1})
        neu_out = true
        muts_AAneu_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        muts_Aaneu_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        muts_aaneu_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(muts_AAneu_next,NaN)
        fill!(muts_Aaneu_next,NaN)
        fill!(muts_aaneu_next,NaN)
    end
    
    for deme in next_gen_posits
        ms1_at_pos = pnt_wld_ms1[deme...]
        ms2_at_pos = pnt_wld_ms2[deme...]

        fitns = []
        cnt_res_AAsel,cnt_res_Aasel,cnt_res_aasel,cnt_res_AAneu,cnt_res_Aaneu,cnt_res_aaneu,fitns =
            calc_muts_and_fitn_in_deme(ms1_at_pos,ms2_at_pos,pnt_wld_stats["s_sel_coef"],pnt_wld_stats["h_domin_coef"],pnt_wld_stats["n_loci"],pnt_wld_stats["sel_loci"])

        if fitn_out
            mean_fitn_next[deme...] = mean(fitns)
        end
        if sel_out
            muts_AAsel_next[deme...] = cnt_res_AAsel
            muts_Aasel_next[deme...] = cnt_res_Aasel
            muts_aasel_next[deme...] = cnt_res_aasel
        end
        if neu_out
            muts_AAneu_next[deme...] = cnt_res_AAneu
            muts_Aaneu_next[deme...] = cnt_res_Aaneu
            muts_aaneu_next[deme...] = cnt_res_aaneu
        end

        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            birth_count = 0
            for _ in 1:next_generation_size
                
                mom_ms1 = wsample(ms1_at_pos,Float32.(fitns))
                mom_ms2 = wsample(ms2_at_pos,Float32.(fitns))
                dad_ms1 = wsample(ms1_at_pos,Float32.(fitns))
                dad_ms2 = wsample(ms2_at_pos,Float32.(fitns))

                gamete_mom_ms1 = copy(mom_ms1)
                gamete_dad_ms1 = copy(dad_ms1)
                gamete_mom_ms2 = copy(mom_ms2)
                gamete_dad_ms2 = copy(dad_ms2)

                crossover(gamete_mom_ms1,gamete_mom_ms2,pnt_wld_stats["n_loci"])
                crossover(gamete_dad_ms1,gamete_dad_ms2,pnt_wld_stats["n_loci"])
                mutate(gamete_mom_ms1,gamete_mom_ms2,pnt_wld_stats["mut_rate"],pnt_wld_stats["n_loci"])
                mutate(gamete_dad_ms1,gamete_dad_ms2,pnt_wld_stats["mut_rate"],pnt_wld_stats["n_loci"])

                move = calc_migr_dist(deme,pnt_wld_stats,migr_mode,bottleneck,max_migr,refl_walls,r_max_migr)

                indices = [deme[i] + move[i] for i in 1:wlddim]

                if !isassigned(wld_ms1_next,indices...)
                    wld_ms1_next[indices...] = []
                    wld_ms2_next[indices...] = []
                end
                push!(wld_ms1_next[indices...],gamete_mom_ms1)
                push!(wld_ms2_next[indices...],gamete_dad_ms2)

                birth_count += 1
                all_birth_count += 1
            end
            
            if pops_out
                pops_next[deme...] = birth_count
            end
        end
    end
    
    #pnt_wld_ms1 = wld_ms1_next
    #pnt_wld_ms2 = wld_ms2_next
    return wld_ms1_next,wld_ms2_next,mean_fitn_next, pops_next, muts_AAsel_next, muts_Aasel_next, muts_aasel_next, muts_AAneu_next, muts_Aaneu_next, muts_aaneu_next
end

@inbounds function create_empty_world(max=(DEF_X_MAX,DEF_Y_MAX);name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),k_capacity=DEF_K_CAPACITY,
    r_prolif_rate=DEF_R_PROLIF_RATE,n_loci=DEF_N_LOCI,n_sel_loci=DEF_N_SEL_LOCI,
    mut_rate=DEF_MUT_RATE,migr_rate=DEF_MIGR_RATE,migr_mode=DEF_MIGR_MODE,s_sel_coef=DEF_S_SEL_COEF,h_domin_coef=DEF_H_DOMIN_COEF,prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

    wld_ms1 = Array{Array{Array{Bool}}}(undef,max...) # array of left (in a pair) monosomes ("ms") of all individuals in space
    wld_ms2 = Array{Array{Array{Bool}}}(undef,max...) # array of right (in a pair) monosomes ("ms") of all individuals in space
    for k in Iterators.product([1:n for n in max]...)
        wld_ms1[k...] = Array{Bool, 1}[]
        wld_ms2[k...] = Array{Bool, 1}[]
    end

    wld_stats = Dict(
        "name" => name,
        "max" => max,
        "k_capacity" => k_capacity,
        "r_prolif_rate" => r_prolif_rate,
        "n_loci" => n_loci,
        "n_sel_loci" => n_sel_loci,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "s_sel_coef" => s_sel_coef,
        "h_domin_coef" => h_domin_coef,
        "prop_of_del_muts" => prop_of_del_muts,
        "wlddim" => length(max)
        #"rangeexps" => []
    )

    return wld_ms1,wld_ms2,wld_stats
end

@inbounds function fill_random_demes(pnt_wld_ms1,pnt_wld_ms2,pnt_wld_stats,max_fill::Vector{UnitRange{Int64}},n_demes_to_fill=DEF_N_DEMES_STARTFILL) # ::Array{Array{Array{Bool}}}
    
    possible_init_coords = [collect(x) for x in Iterators.product(max_fill...)]
    init_coords = sample(possible_init_coords,n_demes_to_fill;replace=false)
    pnt_wld_stats["sel_loci"] = randperm(pnt_wld_stats["n_loci"])[1:pnt_wld_stats["n_sel_loci"]]

    for coord in init_coords
        if !isassigned(pnt_wld_ms1,coord...)
            pnt_wld_ms1[coord...] = []
            pnt_wld_ms2[coord...] = []
        end
        for _ in 1:pnt_wld_stats["k_capacity"]
            push!(pnt_wld_ms1[coord...],falses(pnt_wld_stats["n_loci"]))
            push!(pnt_wld_ms2[coord...],falses(pnt_wld_stats["n_loci"]))
        end
    end

    pnt_wld_stats["startfill"] = copy(max_fill)
    pnt_wld_stats["n_demes_startfill"] = n_demes_to_fill
end

"""
`max` - space bounds in case of creating a new world
"""
function rangeexp(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;max_burnin=(DEF_X_MAX_BURNIN,DEF_Y_MAX),max_exp=(DEF_X_MAX_EXP,DEF_Y_MAX),max=(DEF_X_MAX,DEF_Y_MAX),migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld_ms1=NaN,wld_ms2=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=NaN,r_max_burnin=0,r_max_exp=0)

    fitn_wld = NaN
    pops_wld = NaN
    muts_AAsel_wld = NaN
    muts_Aasel_wld = NaN
    muts_aasel_wld = NaN
    muts_AAneu_wld = NaN
    muts_Aaneu_wld = NaN
    muts_aaneu_wld = NaN    

    if !(wld_ms1 isa Array{Array{Array{Bool}}})
        #println("No world provided. Creating a new world.")
        wld_ms1,wld_ms2,wld_stats = create_empty_world(max;name=name)
        if isa(max_burnin,Tuple)
            startfill_range = [1:upper for upper in max_burnin]
        else
            startfill_range = [(trunc(Int,1+round((upper-1)/2)-r_max_burnin*0.666)):(trunc(Int,1+round((upper-1)/2)+r_max_burnin*0.666)) for upper in max]
        end
        fill_random_demes(wld_ms1,wld_ms2,wld_stats,startfill_range)
    end

    if occursin("F", data_to_generate)
        fitn_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end
    if occursin("P", data_to_generate)
        pops_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end
    if occursin("S", data_to_generate)
        muts_AAsel_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
        muts_Aasel_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
        muts_aasel_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end
    if occursin("N", data_to_generate)
        muts_AAneu_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
        muts_Aaneu_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
        muts_aaneu_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end

    wlddim = wld_stats["wlddim"]

    n_gens_total = n_gens_burnin+n_gens_exp
    @inbounds for g in 1:n_gens_total
        if g<=n_gens_burnin
            max_migr = max_burnin
            r_max_migr = r_max_burnin
        else
            max_migr = max_exp
            r_max_migr = r_max_exp
        end
        wld_ms1,wld_ms2,mean_fitn_next,pops_next,muts_AAsel_next,muts_Aasel_next,muts_aasel_next,muts_AAneu_next,
            muts_Aaneu_next,muts_aaneu_next = build_next_gen(wld_ms1,wld_ms2,wld_stats,fitn_wld,pops_wld,muts_AAsel_wld,muts_Aasel_wld,muts_aasel_wld,muts_AAneu_wld,
            muts_Aaneu_wld,muts_aaneu_wld;
            max_migr=max_migr,migr_mode=migr_mode,bottleneck=bottleneck,r_max_migr=r_max_migr)
        if occursin("F", data_to_generate)
            fitn_wld = cat(fitn_wld, mean_fitn_next, dims=wlddim+1)
        end
        if occursin("P", data_to_generate)
            pops_wld = cat(pops_wld, pops_next, dims=wlddim+1)
        end
        if occursin("S", data_to_generate)
            muts_AAsel_wld = cat(muts_AAsel_wld, muts_AAsel_next, dims=wlddim+1)
            muts_Aasel_wld = cat(muts_Aasel_wld, muts_Aasel_next, dims=wlddim+1)
            muts_aasel_wld = cat(muts_aasel_wld, muts_aasel_next, dims=wlddim+1)
        end
        if occursin("N", data_to_generate)
            muts_AAneu_wld = cat(muts_AAneu_wld, muts_AAneu_next, dims=wlddim+1)
            muts_Aaneu_wld = cat(muts_Aaneu_wld, muts_Aaneu_next, dims=wlddim+1)
            muts_aaneu_wld = cat(muts_aaneu_wld, muts_aaneu_next, dims=wlddim+1)
        end
    end

    #= append!(wld_stats["rangeexps"],Dict(
        "x_max_burnin" => x_max_burnin,
        "y_max_burnin" => DEF_Y_MAX,
        "n_gens_burnin" => n_gens_burnin,
        "n_gens_exp" => n_gens_exp,
        "n_gens" => n_gens_total)) =#
    wld_stats["max_burnin"] = max_burnin
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    wld_stats["n_gens"] = n_gens_total
    
    #= res = []
    function needed_data(symb,arr,res)
        if occursin("P", data_to_generate)
            push!(res,Ref(arr))
        end
    end
    needed_data("P",pops_wld) =#
    return Dict("stats"=>wld_stats,"fitn"=>fitn_wld,"pops"=>pops_wld,"AAsel"=>muts_AAsel_wld,"Aasel"=>muts_Aasel_wld,
        "aasel"=>muts_aasel_wld,"AAneu"=>muts_AAneu_wld,"Aaneu"=>muts_Aaneu_wld,"aaneu"=>muts_aaneu_wld)
end

# Simulation functions (infinite-sites)
# ------------------------------------------------

@inbounds function mutate_inf(person,mut_rate,n_segr_regions,s_sel_coef,prop_of_sel_loci)
    get_mutation_random = rand(Poisson(mut_rate))
    muts_delsel = 0
    muts_bensel = 0
    muts_delneu = 0
    muts_benneu = 0

    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:n_segr_regions*2)

        if rand() < mut_rate
            if rand() < prop_of_sel_loci
                person[pos_alter] *= 1 - s_sel_coef
                muts_delsel += 1
            else
                muts_delneu += 1
            end
        else
            if rand() < prop_of_sel_loci
                person[pos_alter] *= 1 + s_sel_coef
                muts_bensel += 1
            else
                muts_benneu += 1
            end
        end
    end
    
    return muts_delsel, muts_bensel, muts_delneu, muts_benneu
end

@inbounds function crossover_inf(person,n_segr_regions)
    for i in 1:n_segr_regions
        lr = rand(1:2)
        person[i] = lr==1 ? person[i] : person[i+n_segr_regions]
    end
end

@inbounds function mate_inf(person1,person2,n_segr_regions)
    lr1 = rand(1:2)==1 ? (1:n_segr_regions) : ((n_segr_regions+1):(n_segr_regions*2))
    lr2 = rand(1:2)==1 ? (1:n_segr_regions) : ((n_segr_regions+1):(n_segr_regions*2))
    return vcat(person1[lr1],person2[lr2])
end

@inbounds function build_next_gen_inf(pnt_wld::Array{Array{Array{Float32}}},pnt_wld_stats,pnt_fitn_wld=NaN,pnt_pops_wld=NaN,pnt_muts_delsel_wld=NaN,
    pnt_muts_bensel_wld=NaN,pnt_muts_delneu_wld=NaN,pnt_muts_benneu_wld=NaN;
    max_migr=NaN,migr_mode=DEF_MIGR_MODE,bottleneck=NaN,refl_walls=false,r_max_migr=0)

    wlddim = pnt_wld_stats["wlddim"]

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(pnt_wld,pnt_wld_stats)
    
    # Define the habitat (world) and the data arrays in the next generation
    wld_next = Array{Array{Array{Float32}}, wlddim}(undef,pnt_wld_stats["max"]...)
    mean_fitn_next = NaN
    pops_next = NaN
    muts_delsel_next = NaN
    muts_bensel_next = NaN
    muts_delneu_next = NaN
    muts_benneu_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    fitn_out = false
    pops_out = false
    sel_out = false
    neu_out = false

    if pnt_fitn_wld isa Array{Float32, wlddim+1}
        fitn_out = true
        mean_fitn_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(mean_fitn_next,NaN)
    end
    if pnt_pops_wld isa Array{Float32, wlddim+1}
        pops_out = true
        pops_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(pops_next,NaN)
    end
    if (pnt_muts_delsel_wld isa Array{Float32, wlddim+1}) && (pnt_muts_bensel_wld isa Array{Float32, wlddim+1})
        sel_out = true
        muts_delsel_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        muts_bensel_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(muts_delsel_next,NaN)
        fill!(muts_bensel_next,NaN)
    end
    if (pnt_muts_delneu_wld isa Array{Float32, wlddim+1}) && (pnt_muts_benneu_wld isa Array{Float32, wlddim+1})
        neu_out = true
        muts_delneu_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        muts_benneu_next = Array{Float32}(undef,pnt_wld_stats["max"]...)
        fill!(muts_delneu_next,NaN)
        fill!(muts_benneu_next,NaN)
    end
    
    
    for deme in next_gen_posits
        inds_at_pos = pnt_wld[deme...]
        fitns = prod.(inds_at_pos)

        if fitn_out
            mean_fitn_next[deme...] = mean(fitns)
        end

        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            birth_count = 0
            for _ in 1:next_generation_size
                mom = wsample(inds_at_pos,fitns)
                dad = wsample(inds_at_pos,fitns)
                
                gamete_mom = copy(mom)
                gamete_dad = copy(dad)

                crossover_inf(gamete_mom,pnt_wld_stats["n_segr_regions"])
                crossover_inf(gamete_dad,pnt_wld_stats["n_segr_regions"])
                mate_result = mate_inf(gamete_mom,gamete_dad,pnt_wld_stats["n_segr_regions"])

                muts_delsel, muts_bensel, muts_delneu, muts_benneu = mutate_inf(mate_result,pnt_wld_stats["mut_rate"],pnt_wld_stats["n_segr_regions"],pnt_wld_stats["s_sel_coef"],pnt_wld_stats["prop_of_sel_loci"])
                if sel_out
                    if isnan(muts_delsel_next[deme...])
                        muts_delsel_next[deme...] = 0
                    end
                    if isnan(muts_bensel_next[deme...])
                        muts_bensel_next[deme...] = 0
                    end
                    muts_delsel_next[deme...] += muts_delsel
                    muts_bensel_next[deme...] += muts_bensel
                end
                if neu_out
                    if isnan(muts_delneu_next[deme...])
                        muts_delneu_next[deme...] = 0
                    end
                    if isnan(muts_benneu_next[deme...])
                        muts_benneu_next[deme...] = 0
                    end
                    muts_delneu_next[deme...] += muts_delneu
                    muts_benneu_next[deme...] += muts_benneu
                end
                
                move = calc_migr_dist(deme,pnt_wld_stats,migr_mode,bottleneck,max_migr,refl_walls,r_max_migr)

                indices = [deme[i] + move[i] for i in 1:wlddim]

                if !isassigned(wld_next,indices...)
                    wld_next[indices...] = []
                end
                push!(wld_next[indices...],mate_result)

                birth_count += 1
                all_birth_count += 1
            end
            
            if pops_out
                pops_next[deme...] = birth_count
            end
        end
    end
    
    return wld_next, mean_fitn_next, pops_next, muts_delsel_next, muts_bensel_next, muts_delneu_next, muts_benneu_next
end

"""
Create an empty world with infinite-sites individual structure.
"""
function create_empty_world_inf(max=(DEF_X_MAX,DEF_Y_MAX);name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),k_capacity=DEF_K_CAPACITY,
    r_prolif_rate=DEF_R_PROLIF_RATE,n_segr_regions=DEF_N_SEGR_REGIONS,
    mut_rate=DEF_MUT_RATE,migr_rate=DEF_MIGR_RATE,migr_mode=DEF_MIGR_MODE,s_sel_coef=DEF_S_SEL_COEF,prop_of_del_muts=DEF_PROP_OF_DEL_MUTS,prop_of_sel_loci=DEF_PROP_OF_SEL_LOCI)

    wld = Array{Array{Array{Float32}}}(undef,max...) # array of fitness values of all individuals in space
    for k in Iterators.product([1:n for n in max]...)
        wld[k...] = Array{Float32, 1}[]
    end

    wld_stats = Dict(
        "name" => name,
        "max" => max,
        "k_capacity" => k_capacity,
        "r_prolif_rate" => r_prolif_rate,
        "n_segr_regions" => n_segr_regions,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "s_sel_coef" => s_sel_coef,
        "prop_of_del_muts" => prop_of_del_muts,
        "prop_of_sel_loci" => prop_of_sel_loci,
        "wlddim" => length(max)
        #"rangeexps" => []
    )

    return wld,wld_stats
end

function fill_random_demes_inf(pnt_wld::Array{Array{Array{Float32}}},pnt_wld_stats,max_fill::Vector{UnitRange{Int64}},n_demes_to_fill=DEF_N_DEMES_STARTFILL) # ::Array{Array{Array{Bool}}}

    possible_init_coords = [collect(x) for x in Iterators.product(max_fill...)]
    init_coords = sample(possible_init_coords,n_demes_to_fill;replace=false)

    for coord in init_coords
        if !isassigned(pnt_wld,coord...)
            pnt_wld[coord...] = []
        end
        for _ in 1:pnt_wld_stats["k_capacity"]
            push!(pnt_wld[coord...],ones(pnt_wld_stats["n_segr_regions"]*2))
        end
    end

    pnt_wld_stats["startfill"] = copy(max_fill)
    pnt_wld_stats["n_demes_startfill"] = n_demes_to_fill
end

function rangeexp_inf(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;max_burnin=(DEF_X_MAX_BURNIN,DEF_Y_MAX),max_exp=(DEF_X_MAX_EXP,DEF_Y_MAX),max=(DEF_X_MAX,DEF_Y_MAX),migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),prop_of_sel_loci=DEF_PROP_OF_SEL_LOCI,bottleneck=NaN,r_max_burnin=0,r_max_exp=0)

    fitn_wld = NaN
    pops_wld = NaN
    muts_delsel_wld = NaN
    muts_bensel_wld = NaN
    muts_delneu_wld = NaN
    muts_benneu_wld = NaN

    if !(wld isa Array{Array{Array{Float32}}})
        #println("No world provided. Creating a new world.")
        wld,wld_stats = create_empty_world_inf(max;name=name,prop_of_sel_loci=prop_of_sel_loci)
        if isa(max_burnin,Tuple)
            startfill_range = [1:upper for upper in max_burnin]
        else
            startfill_range = [(trunc(Int,1+round((upper-1)/2)-r_max_burnin*0.666)):(trunc(Int,1+round((upper-1)/2)+r_max_burnin*0.666)) for upper in max]
        end
        fill_random_demes_inf(wld,wld_stats,startfill_range)
    end

    if occursin("F", data_to_generate)
        fitn_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end
    if occursin("P", data_to_generate)
        pops_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end
    if occursin("S", data_to_generate)
        muts_delsel_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
        muts_bensel_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end
    if occursin("N", data_to_generate)
        muts_delneu_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
        muts_benneu_wld = Array{Float32}(undef, wld_stats["max"]..., 0)
    end

    wlddim = wld_stats["wlddim"]

    n_gens_total = n_gens_burnin+n_gens_exp
    @inbounds for g in 1:n_gens_total
        if g<=n_gens_burnin
            max_migr = max_burnin
            r_max_migr = r_max_burnin
        else
            max_migr = max_exp
            r_max_migr = r_max_exp
        end

        wld, mean_fitn_next, pops_next, muts_delsel_next, muts_bensel_next, muts_delneu_next, muts_benneu_next = build_next_gen_inf(wld,wld_stats,fitn_wld,pops_wld,muts_delsel_wld,
            muts_bensel_wld,muts_delneu_wld,muts_benneu_wld;
            max_migr=max_migr,migr_mode=migr_mode,bottleneck=bottleneck,r_max_migr=r_max_migr)

        if occursin("F", data_to_generate)
            fitn_wld = cat(fitn_wld, mean_fitn_next, dims=wlddim+1)
        end
        if occursin("P", data_to_generate)
            pops_wld = cat(pops_wld, pops_next, dims=wlddim+1)
        end
        if occursin("S", data_to_generate)
            muts_delsel_wld = cat(muts_delsel_wld, muts_delsel_next, dims=wlddim+1)
            muts_bensel_wld = cat(muts_bensel_wld, muts_bensel_next, dims=wlddim+1)
        end
        if occursin("N", data_to_generate)
            muts_delneu_wld = cat(muts_delneu_wld, muts_delneu_next, dims=wlddim+1)
            muts_benneu_wld = cat(muts_benneu_wld, muts_benneu_next, dims=wlddim+1)
        end
    end

    wld_stats["max_burnin"] = max_burnin
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    wld_stats["n_gens"] = n_gens_total
    
    return Dict("stats"=>wld_stats,"fitn"=>fitn_wld,"pops"=>pops_wld,"delsel"=>muts_delsel_wld,"delneu"=>muts_delneu_wld,
        "bensel"=>muts_bensel_wld,"benneu"=>muts_benneu_wld)
end


# Plotting functions
# ------------------------------------------------

function multi_index(array, fixed_index, fixed_dim)
    dims = ndims(array)  # Number of dimensions of the array
    indices = ntuple(i -> (i == fixed_dim ? fixed_index : :), dims)
    return array[indices...]
end

"""
Shows an animated heatmap of `data` from `gen_start` to `gen_end`.

---

`data`: 3-dimensional array of 2d by-deme data + time 

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`clim`: color bounds (Plots.jl's clim parameter)

`log_base`: if not **-1**, color shows log values with this as base

---

"""
function re_heatmap(re::Dict,dataname::String,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,log_base=-1, kwargs...)
    if !isa(re[dataname],Array)
        println("This data was not generated.")
    else
        re_heatmap(re[dataname],re["stats"]["wlddim"]+1,gen_start,gen_end;slow_factor=slow_factor,log_base=log_base, kwargs...)
    end
end
function re_heatmap(data::Array,dims::Int,gen_start=1,gen_end=DEF_N_GENS_BURNIN+DEF_N_GENS_EXP;slow_factor=1,log_base=-1, kwargs...)
    @gif for i in gen_start:(gen_end*slow_factor-1)
        gen_no = trunc(Int,i/slow_factor)+1

        if all(isnan,multi_index(data,gen_no,dims))
            println("No values found in any deme.")
        end

        if log_base>0 && log_base==1
            obj = log.(log_base,multi_index(data,gen_no,dims)')
        else
            obj = multi_index(data,gen_no,dims)'
        end

        if dims==2 # Including time
            heatmap(obj,ylabel="Generation $gen_no",size=(1200,200); kwargs...)
        else
            heatmap(obj,ylabel="Generation $gen_no"; kwargs...)
        end
        
    end
end


"""
Shows population data of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`clim`: color bounds (Plots.jl's clim parameter)

`log_base`: if not **-1**, color shows log values with this as base

---

"""
function re_heatmap_pops(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["k_capacity"]),log_base=-1,plot_options...)
    re_heatmap(re,"pops",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_fitn(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,1),log_base=-1,plot_options...)
    re_heatmap(re,"fitn",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_AAsel(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,length(re["stats"]["sel_loci"])),log_base=-1,plot_options...)
    re_heatmap(re,"AAsel",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_Aasel(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,length(re["stats"]["sel_loci"])),log_base=-1,plot_options...)
    re_heatmap(re,"Aasel",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_aasel(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,length(re["stats"]["sel_loci"])),log_base=-1,plot_options...)
    re_heatmap(re,"aasel",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_AAneu(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1,plot_options...)
    re_heatmap(re,"AAneu",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_Aaneu(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1,plot_options...)
    re_heatmap(re,"Aaneu",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_aaneu(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1,plot_options...)
    re_heatmap(re,"aaneu",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_delsel(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1,plot_options...)
    re_heatmap(re,"delsel",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_delneu(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1,plot_options...)
    re_heatmap(re,"delneu",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_bensel(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1,plot_options...)
    re_heatmap(re,"bensel",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

function re_heatmap_benneu(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1,plot_options...)
    re_heatmap(re,"benneu",gen_start,gen_end;slow_factor=slow_factor,clim=clim,log_base=log_base,plot_options...)
end

# Functions pertaining to averaging and expansion front
# ------------------------------------------------

function average_all(data::Array,n_gens::Int;greaterzero=false,divide=true)
    res = Array{Float32}(undef,0)
    for j in 1:n_gens
        push!(res,mean(skipmissing(multi_index(data,j,dims))))
    end
    return res
end
function average_all(re::Dict,dataname::String;greaterzero=false,divide=true)
    res = Array{Float32}(undef,0)
    for j in 1:re["stats"]["n_gens"]
        push!(res,mean(filter(!isnan,multi_index(re[dataname],j,dims))))
    end
    return res
end

"""
Finds the average value of `dataname` between all demes at the expansion front of `re`.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`greaterzero`: if **true**, **>0** values are considered when determining the front (**>=0** values if **false**)

`oneside`: if **true**, approach only from one side (i.e. from the positive direction in strip expansions)

`divide`: if **true**, find average

---

"""
function average_front(re,dataname;greaterzero=false,oneside=false,divide=true)
    average_front(re[dataname],re["stats"]["n_gens"],re["stats"]["max"]...;greaterzero=greaterzero,oneside=oneside,divide=divide)
end

function average_front(data::Array,n_gens,x_max,y_max;greaterzero=false,oneside=false,divide=true)
    front_array = Array{Float32}(undef,0)
    for j in 1:n_gens
        sum_total = 0
        cnt = 0
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && (isnan(data[frontier_x,_y,j]) || (greaterzero && data[frontier_x,_y,j] == 0))
                frontier_x -= 1
            end
            if data[frontier_x,_y,j]>=0 || (greaterzero && data[frontier_x,_y,j]>0)
                sum_total += data[frontier_x,_y,j]
                cnt += 1
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && (isnan(data[frontier_x,_y,j]) || (greaterzero && data[frontier_x,_y,j] == 0))
                    frontier_x += 1
                end
                if data[frontier_x,_y,j]>=0 || (greaterzero && data[frontier_x,_y,j]>0)
                    sum_total += data[frontier_x,_y,j]
                    cnt += 1
                end
            end
        end
        mean_both_sides_y = sum_total
        if divide
            mean_both_sides_y /= cnt
        end

        if !oneside
            sum_total = 0
            cnt = 0
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && (isnan(data[_x,frontier_y,j]) || (greaterzero && data[_x,frontier_y,j] == 0))
                    frontier_y -= 1
                end
                if data[_x,frontier_y,j]>=0 || (greaterzero && data[_x,frontier_y,j]>0)
                    sum_total += data[_x,frontier_y,j]
                    cnt += 1
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && (isnan(data[_x,frontier_y,j]) || (greaterzero && data[_x,frontier_y,j] == 0))
                    frontier_y += 1
                end
                if data[_x,frontier_y,j]>=0 || (greaterzero && data[_x,frontier_y,j]>0)
                    sum_total += data[_x,frontier_y,j]
                    cnt += 1
                end
            end
            if divide
                mean_both_sides_x = sum_total/cnt
            end
            front_array = cat(front_array,(mean_both_sides_x+mean_both_sides_y)/2, dims=1)
        else
            front_array = cat(front_array,mean_both_sides_y, dims=1)
        end
    end
    return front_array
end

"""
Finds the front array of `dataname` in `re`.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`oneside`: if **true**, approach only from one side (i.e. from the positive direction in strip expansions)

---

"""
function front_array(re,dataname;oneside=false)
    front_array(re[dataname],re["stats"]["n_gens"],re["stats"]["max"]...;oneside=oneside)
end

# 1D
function front_array(data,n_gens,x_max;oneside=false)
    front_arr = fill(NaN,x_max,n_gens)

    for j in 1:n_gens
        frontier = x_max
        while frontier != 1 && isnan(data[frontier,j])
            frontier -= 1
        end
        if !isnan(data[frontier,j])
            front_arr[frontier,j]=data[frontier,j]
        end
        if !oneside
            frontier = 1
            while frontier != x_max && isnan(data[frontier,j])
                frontier += 1
            end
            if !isnan(data[frontier,j])
                front_arr[frontier,j]=data[frontier,j]
            end
        end
    end
    return front_arr
end

# 2D
function front_array(data::Array,n_gens,x_max,y_max;oneside=false)
    front_arr = fill(NaN,x_max,y_max,n_gens)
    for j in 1:n_gens
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && isnan(data[frontier_x,_y,j])
                frontier_x -= 1
            end
            if !isnan(data[frontier_x,_y,j])
                front_arr[frontier_x,_y,j]=data[frontier_x,_y,j]
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && isnan(data[frontier_x,_y,j])
                    frontier_x += 1
                end
                if !isnan(data[frontier_x,_y,j])
                    front_arr[frontier_x,_y,j]=data[frontier_x,_y,j]
                end
            end
        end

        if !oneside
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && isnan(data[_x,frontier_y,j])
                    frontier_y -= 1
                end
                if !isnan(data[_x,frontier_y,j])
                    front_arr[_x,frontier_y,j]=data[_x,frontier_y,j]
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && isnan(data[_x,frontier_y,j])
                    frontier_y += 1
                end
                if !isnan(data[_x,frontier_y,j]>0)
                    front_arr[_x,frontier_y,j]=data[_x,frontier_y,j]
                end
            end
        end
    end
    return front_arr
end

# mean front fitness (or other data)
"""
Normalises `dataname` in `re` using the "maximum normalisation" method: after the last burn-in generation, divide by the maximum of each generation.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

---

"""
function norm_maximum(re,dataname)
    norm_maximum(re[dataname],re["stats"]["wlddim"]+1,re["stats"]["n_gens_burnin"],re["stats"]["n_gens_exp"])
end

function norm_maximum(data::Array,dims::Int,n_gens_burnin,n_gens_exp)
    normal_array = copy(data)
    for j in 1:n_gens_exp
        gen_max = maximum(multi_index(data,n_gens_burnin+j,dims))
        multi_index(normal_array,n_gens_burnin+j,dims) ./= gen_max
    end
    return normal_array
end

"""
Normalises `dataname` in `re` using the "maximum normalisation" method: after the last burn-in generation, divide by the constant value of average fitness over all demes at the (last burn-in generation+1)=onset generation

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`offset` - offset from the onset generation

---

"""
function norm_onset_mean(re::Dict,dataname::String="fitn",offset=0)
    norm_onset_mean(re[dataname],re["stats"]["wlddim"]+1,re["stats"]["n_gens_burnin"],offset)
end

function norm_onset_mean(data::Array,dims::Int,n_gens_burnin::Int,offset=0)
    normal_array = copy(data)

    sum = 0
    count = 0
    for u in multi_index(data,n_gens_burnin+1+offset,dims)
        if u > 0
            sum += u
            count += 1
        end
    end
    gen_average = sum/count

    for i in (n_gens_burnin+1):length(normal_array)
        multi_index(data,i,dims) ./= gen_average
    end
    return normal_array
end

#= function normalise_front_by_onset_mean(average_1d_array)
    normal_array = copy(average_1d_array)
    normal_array[BURN_IN_GEN_N+1:end] /= average_1d_array[BURN_IN_GEN_N+1]
    return normal_array
end

function normalise_front_by_max(average_1d_array,fitn_array)
    normal_array = copy(average_1d_array)
    for j in 1:(TOTAL_GEN_N-BURN_IN_GEN_N)
        gen_max = maximum(fitn_array[:,:,BURN_IN_GEN_N+j])
        normal_array[BURN_IN_GEN_N+j] /= gen_max
    end
    return normal_array
end

function find_front_array_muts(data,muts_array;oneside=false)
    res_muts = zeros(Float32,y_max,n_gen)
    for j in 1:n_gen
        # scanning every y: side 1

        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && data[frontier_x,_y,j] < 0
                frontier_x -= 1
            end
            if data[frontier_x,_y,j]>0
                res_muts[_y,j]=muts_array[frontier_x,_y,j]
            end
        end
        # scanning every y: side 2
        # add later
    end
    return res_muts
end =#




# 1D
# ------------------------------------------------

function rangeexp_1d(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld_ms1=NaN,wld_ms2=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=NaN)

    rangeexp(n_gens_burnin,n_gens_exp;max_burnin=(x_max_burnin,),max_exp=(x_max_exp,),max=(x_max_exp,),
        migr_mode=migr_mode,data_to_generate=data_to_generate,wld_ms1=wld_ms1,wld_ms2=wld_ms2,wld_stats=wld_stats,name=name,bottleneck=bottleneck)
end


# 2D
# ------------------------------------------------

"""
Simulates an strip range expansion, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in

`x_max_exp`: the outward x-coordinate bound for migration during the expansion

`y_max`: the upper y-coordinate bound (lower bound is always **1** currently)

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**)
- **M** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**)

If starting from existing world, also provide:

`wld_ms1`: world left monosome array

`wld_ms2`: world right monosome array

`wld_stats`: world stats Dict

---
"""
function rangeexp_strip(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,y_max=DEF_Y_MAX,migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld_ms1=NaN,wld_ms2=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=("midhole at x=",x_max_burnin*2))

    rangeexp(n_gens_burnin,n_gens_exp;max_burnin=(x_max_burnin,y_max),max_exp=(x_max_exp,y_max),max=(x_max_exp,y_max),
        migr_mode=migr_mode,data_to_generate=data_to_generate,wld_ms1=wld_ms1,wld_ms2=wld_ms2,wld_stats=wld_stats,name=name,bottleneck=bottleneck)
end


"""
Simulates an strip range expansion, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in

`x_max_exp`: the outward x-coordinate bound for migration during the expansion

`y_max`: the upper y-coordinate bound (lower bound is always **0** currently)

`migr_mode`: mode of migration. Possible values:
- **4** - orthogonal directions only
- **6** - hexagonal grid
- **8** - orthogonal and diagonal
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of selected mutations (**sel**)
- **M** - deme-average number of neutral mutations (**neu**)

If starting from existing world, also provide:

`wld`: world fitness array

`wld_stats`: world stats Dict

---

"""
function rangeexp_strip_inf(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,y_max=DEF_Y_MAX,migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=("midhole at x=",x_max_burnin*2),prop_of_sel_loci=DEF_PROP_OF_SEL_LOCI)

    rangeexp_inf(n_gens_burnin,n_gens_exp;max_burnin=(x_max_burnin,y_max),max_exp=(x_max_exp,y_max),max=(x_max_exp,y_max),
        migr_mode=migr_mode,data_to_generate=data_to_generate,wld=wld,wld_stats=wld_stats,name=name,prop_of_sel_loci=prop_of_sel_loci,bottleneck=bottleneck)
end

function rangeexp_disk(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;r_max_burnin=DEF_R_MAX_BURNIN,r_max_exp=DEF_R_MAX_EXP,migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld_ms1=NaN,wld_ms2=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=NaN)

    rangeexp(n_gens_burnin,n_gens_exp;r_max_burnin=r_max_burnin,r_max_exp=r_max_exp,max_burnin=NaN,max_exp=NaN,max=(r_max_exp*2+1,r_max_exp*2+1),
        migr_mode=migr_mode,data_to_generate=data_to_generate,wld_ms1=wld_ms1,wld_ms2=wld_ms2,wld_stats=wld_stats,name=name,bottleneck=bottleneck)
end

function rangeexp_disk_inf(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;r_max_burnin=DEF_R_MAX_BURNIN,r_max_exp=DEF_R_MAX_EXP,migr_mode=DEF_MIGR_MODE,max=(r_max_exp*2+1,r_max_exp*2+1),
    data_to_generate=DEF_DATA_TO_GENERATE,wld=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=NaN,max_exp=NaN)

    rangeexp_inf(n_gens_burnin,n_gens_exp;r_max_burnin=r_max_burnin,r_max_exp=r_max_exp,max_burnin=NaN,max_exp=max_exp,max=max,
        migr_mode=migr_mode,data_to_generate=data_to_generate,wld=wld,wld_stats=wld_stats,name=name,bottleneck=bottleneck)
end

# Future work

vc(x) = cat(eachslice(x, dims=4)...,dims=2)

function re_get_avrel(data::Array,x,gen,denom)
    nd = ndims(data)
    if nd==4
        return mean(vc(data)[x,:,gen])/denom
    elseif nd==3
        return mean(data[x,:,gen])/denom
    else
        println("Wrong data type.")
    end
end
function re_get_avrel(re::Dict,dataname::String,x,gen=Int(re["stats"]["n_gens"]);sel=true)
    denom = sel ? re["stats"]["n_sel_loci"] : re["stats"]["n_loci"]-re["stats"]["n_sel_loci"]
    return re_get_avrel(re[dataname],x,gen,denom)
end
function re_get_avrelAAsel(re::Dict,x,gen=re["stats"]["n_gens"])
    return re_get_avrel(re,"AAsel",x,gen;sel=true)
end
function re_get_avrelAasel(re::Dict,x,gen=re["stats"]["n_gens"])
    return re_get_avrel(re,"Aasel",x,gen;sel=true)
end
function re_get_avrelaasel(re::Dict,x,gen=re["stats"]["n_gens"])
    return re_get_avrel(re,"aasel",x,gen;sel=true)
end
function re_get_avrelAAneu(re::Dict,x,gen=re["stats"]["n_gens"])
    return re_get_avrel(re,"AAneu",x,gen;sel=false)
end
function re_get_avrelAaneu(re::Dict,x,gen=re["stats"]["n_gens"])
    return re_get_avrel(re,"Aaneu",x,gen;sel=false)
end
function re_get_avrelaaneu(re::Dict,x,gen=re["stats"]["n_gens"])
    return re_get_avrel(re,"aaneu",x,gen;sel=false)
end

function re_plot_avrelselneu(re::Dict,dataname::String,x_range=(1:Int(re["stats"]["x_max"]));x_scale_factor=1,sel=true,overlay=false)
    nd = ndims(re[dataname*"sel"])
    if nd==4
        data1 = vc(re[dataname*"sel"])
        data2 = vc(re[dataname*"neu"])
    else
        data1 = re[dataname*"sel"]
        data2 = re[dataname*"neu"]
    end
    t = [re_get_avrel(data1,j,Int(re["stats"]["n_gens"]),re["stats"]["n_sel_loci"]) for j in x_range]
    t2 = [re_get_avrel(data2,j,Int(re["stats"]["n_gens"]),re["stats"]["n_loci"]-re["stats"]["n_sel_loci"]) for j in x_range]

    if haskey(re["stats"],"name")
        lbl1 = re["stats"]["name"]*"[selected $dataname]"
        lbl2 = re["stats"]["name"]*"[neutral $dataname]"
    else
        lbl1 = "selected $dataname"
        lbl2 = "neutral $dataname"
    end

    if overlay
        plot!(x_range*x_scale_factor,t,label=lbl1,xlabel="x")
    else
        plot(x_range*x_scale_factor,t,label=lbl1,xlabel="x")
    end
    #plot!(x_range*x_scale_factor,t2,label=lbl2)
end

function re_plot_avrelselneu!(re::Dict,dataname::String,x_range=(1:Int(re["stats"]["x_max"]));x_scale_factor=1,sel=true,overlay=false)
    re_plot_avrelselneu(re,dataname,x_range;x_scale_factor=x_scale_factor,sel=sel,overlay=true)
end

# 3D
# ------------------------------------------------

function rangeexp_3d_cyl(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld_ms1=NaN,wld_ms2=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),bottleneck=NaN)

    rangeexp(n_gens_burnin,n_gens_exp;max_burnin=(x_max_burnin,),max_exp=(x_max_exp,),max=(x_max_exp,),
        migr_mode=migr_mode,data_to_generate=data_to_generate,wld_ms1=wld_ms1,wld_ms2=wld_ms2,wld_stats=wld_stats,name=name,bottleneck=bottleneck)
end




println("RESK successfully loaded.")