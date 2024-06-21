include("core.jl")

# Simulation functions (regular)
# ------------------------------------------------

@inbounds function calc_migr_dist_1d(deme,pnt_wld_stats,migr_mode,x_bottleneck,x_max_migr,refl_walls)

    move_x = 0
    migr_res = rand()
    if rand() < pnt_wld_stats["migr_rate"] && migr_res < p_lat + p_diag

        # Raw migration results
        move_x = dir[1]
        move_y = dir[2]

        # Nullify migration on certain conditions
        if !isnan(x_bottleneck) && (deme[1]+move_x==x_bottleneck)
            move_x = 0
            move_y = 0
        else
            if deme[1]+move_x > x_max_migr || deme[1]+move_x < 1 # burn-in area check
                move_x = refl_walls ? -move_x : 0
            end
        end
    end
    return move_x
end

@inbounds function build_next_gen(pnt_wld_ms1::Array{Array{Array{Bool}}},pnt_wld_ms2::Array{Array{Array{Bool}}},pnt_wld_stats,pnt_meanf_wld=NaN,pnt_pops_wld=NaN,pnt_muts_AAsel_wld=NaN,pnt_muts_Aasel_wld=NaN,
    pnt_muts_aasel_wld=NaN,pnt_muts_AAneu_wld=NaN,pnt_muts_Aaneu_wld=NaN,pnt_muts_aaneu_wld=NaN;
    x_max_migr=NaN,y_max_migr=NaN,migr_mode=DEF_MIGR_MODE_2D,x_bottleneck=NaN,refl_walls=false)

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(pnt_wld_ms1,pnt_wld_stats)
    
    # Define the habitat (world) and the data arrays in the next generation
    wld_ms1_next = Array{Array{Array{Bool}}}(undef,pnt_wld_stats["x_max"])
    wld_ms2_next = Array{Array{Array{Bool}}}(undef,pnt_wld_stats["x_max"])
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
    meanf_out = false
    pops_out = false
    sel_out = false
    neu_out = false

    if pnt_meanf_wld isa Array{Float32, 3}
        meanf_out = true
        mean_fitn_next = Array{Float32}(undef,pnt_wld_stats["x_max"])
        fill!(mean_fitn_next,NaN)
    end
    if pnt_pops_wld isa Array{Float32, 3}
        pops_out = true
        pops_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(pops_next,NaN)
    end
    if (pnt_muts_AAsel_wld isa Array{Float32, 3}) && (pnt_muts_Aasel_wld isa Array{Float32, 3}) && (pnt_muts_aasel_wld isa Array{Float32, 3})
        sel_out = true
        muts_AAsel_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        muts_Aasel_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        muts_aasel_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(muts_AAsel_next,NaN)
        fill!(muts_Aasel_next,NaN)
        fill!(muts_aasel_next,NaN)
    end
    if (pnt_muts_AAneu_wld isa Array{Float32, 3}) && (pnt_muts_Aaneu_wld isa Array{Float32, 3}) && (pnt_muts_aaneu_wld isa Array{Float32, 3})
        neu_out = true
        muts_AAneu_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        muts_Aaneu_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        muts_aaneu_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(muts_AAneu_next,NaN)
        fill!(muts_Aaneu_next,NaN)
        fill!(muts_aaneu_next,NaN)
    end
    
    for deme in next_gen_posits
        ms1_at_pos = pnt_wld_ms1[deme...]
        ms2_at_pos = pnt_wld_ms2[deme...]

        fitns = []
        cnt_res_AAsel,cnt_res_Aasel,cnt_res_aasel,cnt_res_AAneu,cnt_res_Aaneu,cnt_res_aaneu,fitns =
            calc_muts_and_meanf_in_deme(ms1_at_pos,ms2_at_pos,pnt_wld_stats["s_sel_coef"],pnt_wld_stats["h_domin_coef"],pnt_wld_stats["n_loci"],pnt_wld_stats["sel_loci"])

        if meanf_out
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

                move_x, move_y = calc_migr_dist(deme,pnt_wld_stats,migr_mode,x_bottleneck,x_max_migr,refl_walls)

                if !isassigned(wld_ms1_next,deme[1]+move_x,deme[2]+move_y)
                    wld_ms1_next[deme[1]+move_x,deme[2]+move_y] = []
                    wld_ms2_next[deme[1]+move_x,deme[2]+move_y] = []
                end
                push!(wld_ms1_next[deme[1]+move_x,deme[2]+move_y],gamete_mom_ms1)
                push!(wld_ms2_next[deme[1]+move_x,deme[2]+move_y],gamete_dad_ms2)

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

@inbounds function create_empty_world(x_max=DEF_X_MAX,y_max=DEF_Y_MAX;name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),k_capacity=DEF_K_CAPACITY,
    r_prolif_rate=DEF_R_PROLIF_RATE,n_loci=DEF_N_LOCI,n_sel_loci=DEF_N_SEL_LOCI,
    mut_rate=DEF_MUT_RATE,migr_rate=DEF_MIGR_RATE,migr_mode=DEF_MIGR_MODE_2D,s_sel_coef=DEF_S_SEL_COEF,h_domin_coef=DEF_H_DOMIN_COEF,prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

    wld_ms1 = Array{Array{Array{Bool}}}(undef,x_max,y_max) # array of left (in a pair) monosomes ("ms") of all individuals in space
    wld_ms2 = Array{Array{Array{Bool}}}(undef,x_max,y_max) # array of right (in a pair) monosomes ("ms") of all individuals in space
    for i in 1:x_max, j in 1:y_max
        wld_ms1[i, j] = Array{Bool, 1}[]
        wld_ms2[i, j] = Array{Bool, 1}[]
    end

    wld_stats = Dict(
        "name" => name,
        "x_max" => x_max,
        "y_max" => y_max,
        "k_capacity" => k_capacity,
        "r_prolif_rate" => r_prolif_rate,
        "n_loci" => n_loci,
        "n_sel_loci" => n_sel_loci,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "s_sel_coef" => s_sel_coef,
        "h_domin_coef" => h_domin_coef,
        "prop_of_del_muts" => prop_of_del_muts

        #"rangeexps" => []
    )

    return wld_ms1,wld_ms2,wld_stats
end

function fill_random_demes(pnt_wld_ms1::Array{Array{Array{Bool}}},pnt_wld_ms2::Array{Array{Array{Bool}}},pnt_wld_stats,x_max_fill,y_max_fill,n_demes_to_fill=DEF_N_DEMES_STARTFILL)

    possible_init_coords = [collect(x) for x in Iterators.product(1:x_max_fill, 1:y_max_fill)]
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

    pnt_wld_stats["x_startfill"] = x_max_fill
    pnt_wld_stats["y_startfill"] = y_max_fill
    pnt_wld_stats["n_demes_startfill"] = n_demes_to_fill
end

"""
Simulates an axial range expansion, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with ```DEF_N_DEMES_STARTFILL``` demes filled with individuals.

---

```n_gens_burnin```: duration of the burn-in phase, used to reach mutation-selection equilibrium

```n_gens_exp```: duration of the expansion

```x_max_burnin```: the outward x-coordinate bound for migration during burn-in

```x_max_exp```: the outward x-coordinate bound for migration during the expansion

```y_max```: the upper y-coordinate bound (lower bound is always **0** currently)

```migr_mode```: mode of migration. Possible values:
- **4** - lateral directions only
- **6** - hexagonal grid
- **8** - lateral and diagonal
- **diag1/2** - lateral and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**meanf**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**)
- **M** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**)

If starting from existing world, also provide:

```wld_ms1```: world left monosome array

```wld_ms2```: world right monosome array

```wld_stats```: world stats Dict

---

"""
function rangeexp_axial(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,y_max=DEF_Y_MAX,migr_mode=DEF_MIGR_MODE_2D,
    data_to_generate=DEF_DATA_TO_GENERATE,wld_ms1=NaN,wld_ms2=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS")) # expansion along x model 2

    meanf_wld = NaN
    pops_wld = NaN
    muts_AAsel_wld = NaN
    muts_Aasel_wld = NaN
    muts_aasel_wld = NaN
    muts_AAneu_wld = NaN
    muts_Aaneu_wld = NaN
    muts_aaneu_wld = NaN

    if !(wld_ms1 isa Array{Array{Array{Bool}}})
        #println("No world provided. Creating a new world.")
        wld_ms1,wld_ms2,wld_stats = create_empty_world(x_max_exp,y_max;name=name)
        fill_random_demes(wld_ms1,wld_ms2,wld_stats,Int(x_max_exp/20),y_max)
    end

    if occursin("F", data_to_generate)
        meanf_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("P", data_to_generate)
        pops_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("S", data_to_generate)
        muts_AAsel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_Aasel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_aasel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("N", data_to_generate)
        muts_AAneu_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_Aaneu_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_aaneu_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end

    n_gens_total = n_gens_burnin+n_gens_exp
    @inbounds for g in 1:n_gens_total
        if g<=n_gens_burnin
            _x_max_used = x_max_burnin
        else
            _x_max_used = x_max_exp
        end
        wld_ms1,wld_ms2,mean_fitn_next,pops_next,muts_AAsel_next,muts_Aasel_next,muts_aasel_next,muts_AAneu_next,
        muts_Aaneu_next,muts_aaneu_next = build_next_gen(wld_ms1,wld_ms2,wld_stats,meanf_wld,pops_wld,muts_AAsel_wld,muts_Aasel_wld,muts_aasel_wld,muts_AAneu_wld,
        muts_Aaneu_wld,muts_aaneu_wld;x_max_migr=_x_max_used,y_max_migr=DEF_Y_MAX,migr_mode=migr_mode,x_bottleneck=x_max_burnin*2)
        if occursin("F", data_to_generate)
            meanf_wld = cat(meanf_wld, mean_fitn_next, dims=3)
        end
        if occursin("P", data_to_generate)
            pops_wld = cat(pops_wld, pops_next, dims=3)
        end
        if occursin("S", data_to_generate)
            muts_AAsel_wld = cat(muts_AAsel_wld, muts_AAsel_next, dims=3)
            muts_Aasel_wld = cat(muts_Aasel_wld, muts_Aasel_next, dims=3)
            muts_aasel_wld = cat(muts_aasel_wld, muts_aasel_next, dims=3)
        end
        if occursin("N", data_to_generate)
            muts_AAneu_wld = cat(muts_AAneu_wld, muts_AAneu_next, dims=3)
            muts_Aaneu_wld = cat(muts_Aaneu_wld, muts_Aaneu_next, dims=3)
            muts_aaneu_wld = cat(muts_aaneu_wld, muts_aaneu_next, dims=3)
        end
    end

    #= append!(wld_stats["rangeexps"],Dict(
        "x_max_burnin" => x_max_burnin,
        "y_max_burnin" => DEF_Y_MAX,
        "n_gens_burnin" => n_gens_burnin,
        "n_gens_exp" => n_gens_exp,
        "n_gens" => n_gens_total)) =#
    wld_stats["x_max_burnin"] = x_max_burnin
    wld_stats["y_max_burnin"] = DEF_Y_MAX
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
    return Dict("stats"=>wld_stats,"meanf"=>meanf_wld,"pops"=>pops_wld,"AAsel"=>muts_AAsel_wld,"Aasel"=>muts_Aasel_wld,
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

@inbounds function build_next_gen_inf(pnt_wld::Array{Array{Array{Float32}}},pnt_wld_stats,pnt_meanf_wld=NaN,pnt_pops_wld=NaN,pnt_muts_delsel_wld=NaN,pnt_muts_bensel_wld=NaN,pnt_muts_delneu_wld=NaN,pnt_muts_benneu_wld=NaN;
    x_max_migr=NaN,y_max_migr=NaN,migr_mode=DEF_MIGR_MODE_2D,x_bottleneck=NaN,refl_walls=false)

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(pnt_wld,pnt_wld_stats)
    
    # Define the habitat (world) and the data arrays in the next generation
    wld_next = Array{Array{Array{Float32}}}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
    mean_fitn_next = NaN
    pops_next = NaN
    muts_delsel_next = NaN
    muts_bensel_next = NaN
    muts_delneu_next = NaN
    muts_benneu_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    meanf_out = false
    pops_out = false
    sel_out = false
    neu_out = false
    if pnt_meanf_wld isa Array{Float32, 3}
        meanf_out = true
        mean_fitn_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(mean_fitn_next,NaN)
    end
    if pnt_pops_wld isa Array{Float32, 3}
        pops_out = true
        pops_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(pops_next,NaN)
    end
    if (pnt_muts_delsel_wld isa Array{Float32, 3}) && (pnt_muts_bensel_wld isa Array{Float32, 3})
        sel_out = true
        muts_delsel_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        muts_bensel_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(muts_delsel_next,NaN)
        fill!(muts_bensel_next,NaN)
    end
    if (pnt_muts_delneu_wld isa Array{Float32, 3}) && (pnt_muts_benneu_wld isa Array{Float32, 3})
        neu_out = true
        muts_delneu_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        muts_benneu_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(muts_delneu_next,NaN)
        fill!(muts_benneu_next,NaN)
    end
    
    
    for deme in next_gen_posits
        inds_at_pos = pnt_wld[deme...]
        fitns = prod.(inds_at_pos)

        if meanf_out
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
                
                move_x, move_y = calc_migr_dist(deme,pnt_wld_stats,migr_mode,x_bottleneck,x_max_migr,refl_walls)

                if !isassigned(wld_next,deme[1]+move_x,deme[2]+move_y)
                    wld_next[deme[1]+move_x,deme[2]+move_y] = []
                end
                push!(wld_next[deme[1]+move_x,deme[2]+move_y],mate_result)

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
function create_empty_world_inf(x_max=DEF_X_MAX,y_max=DEF_Y_MAX;name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),k_capacity=DEF_K_CAPACITY,
    r_prolif_rate=DEF_R_PROLIF_RATE,n_segr_regions=DEF_N_SEGR_REGIONS,n_sel_loci=DEF_N_SEL_LOCI,
    mut_rate=DEF_MUT_RATE,migr_rate=DEF_MIGR_RATE,migr_mode=DEF_MIGR_MODE_2D,s_sel_coef=DEF_S_SEL_COEF,prop_of_del_muts=DEF_PROP_OF_DEL_MUTS,prop_of_sel_loci=DEF_PROP_OF_SEL_LOCI)

    wld = Array{Array{Array{Float32}}}(undef,x_max,y_max) # array of fitness values of all individuals in space
    for i in 1:x_max, j in 1:y_max
        wld[i, j] = Array{Float32, 1}[]
    end

    wld_stats = Dict(
        "name" => name,
        "x_max" => x_max,
        "y_max" => y_max,
        "k_capacity" => k_capacity,
        "r_prolif_rate" => r_prolif_rate,
        "n_segr_regions" => n_segr_regions,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "s_sel_coef" => s_sel_coef,
        "prop_of_del_muts" => prop_of_del_muts,
        "prop_of_sel_loci" => prop_of_sel_loci

        #"rangeexps" => []
    )

    return wld,wld_stats
end

function fill_random_demes_inf(pnt_wld::Array{Array{Array{Float32}}},pnt_wld_stats,x_max_fill,y_max_fill,n_demes_to_fill=DEF_N_DEMES_STARTFILL)

    possible_init_coords = [collect(x) for x in Iterators.product(1:x_max_fill, 1:y_max_fill)]
    init_coords = sample(possible_init_coords,n_demes_to_fill;replace=false)

    for coord in init_coords
        if !isassigned(pnt_wld,coord...)
            pnt_wld[coord...] = []
        end
        for _ in 1:pnt_wld_stats["k_capacity"]
            push!(pnt_wld[coord...],ones(pnt_wld_stats["n_segr_regions"]*2))
        end
    end

    pnt_wld_stats["x_startfill"] = x_max_fill
    pnt_wld_stats["y_startfill"] = y_max_fill
    pnt_wld_stats["n_demes_startfill"] = n_demes_to_fill
end

"""
Simulates an axial range expansion, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with ```DEF_N_DEMES_STARTFILL``` demes filled with individuals.

---

```n_gens_burnin```: duration of the burn-in phase, used to reach mutation-selection equilibrium

```n_gens_exp```: duration of the expansion

```x_max_burnin```: the outward x-coordinate bound for migration during burn-in

```x_max_exp```: the outward x-coordinate bound for migration during the expansion

```y_max```: the upper y-coordinate bound (lower bound is always **0** currently)

```migr_mode```: mode of migration. Possible values:
- **4** - lateral directions only
- **6** - hexagonal grid
- **8** - lateral and diagonal
- **diag1/2** - lateral and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**meanf**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of selected mutations (**sel**)
- **M** - deme-average number of neutral mutations (**neu**)

If starting from existing world, also provide:

```wld```: world fitness array

```wld_stats```: world stats Dict

---

"""
function rangeexp_axial_inf(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,y_max=DEF_Y_MAX,migr_mode=DEF_MIGR_MODE_2D,
    data_to_generate=DEF_DATA_TO_GENERATE,wld=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),prop_of_sel_loci=DEF_PROP_OF_SEL_LOCI) # expansion along x model 2

    meanf_wld = NaN
    pops_wld = NaN
    muts_delsel_wld = NaN
    muts_bensel_wld = NaN
    muts_delneu_wld = NaN
    muts_benneu_wld = NaN

    if !(wld isa Array{Array{Array{Float32}}})
        #println("No world provided. Creating a new world.")
        wld,wld_stats = create_empty_world_inf(x_max_exp,y_max;name=name,prop_of_sel_loci=prop_of_sel_loci)
        fill_random_demes_inf(wld,wld_stats,Int(x_max_exp/20),y_max)
    end

    if occursin("F", data_to_generate)
        meanf_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("P", data_to_generate)
        pops_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("S", data_to_generate)
        muts_delsel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_bensel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("N", data_to_generate)
        muts_delneu_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_benneu_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end

    n_gens_total = n_gens_burnin+n_gens_exp

    @inbounds for g in 1:n_gens_total
        if g<=n_gens_burnin
            _x_max_used = x_max_burnin
        else
            _x_max_used = x_max_exp
        end

        wld, mean_fitn_next, pops_next, muts_delsel_next, muts_bensel_next, muts_delneu_next, muts_benneu_next = build_next_gen_inf(wld,wld_stats,meanf_wld,pops_wld,muts_delsel_wld,
            muts_bensel_wld,muts_delneu_wld,muts_benneu_wld;
            x_max_migr=_x_max_used,y_max_migr=DEF_Y_MAX,migr_mode=migr_mode,x_bottleneck=x_max_burnin*2)

        if occursin("F", data_to_generate)
            meanf_wld = cat(meanf_wld, mean_fitn_next, dims=3)
        end
        if occursin("P", data_to_generate)
            pops_wld = cat(pops_wld, pops_next, dims=3)
        end
        if occursin("S", data_to_generate)
            muts_delsel_wld = cat(muts_delsel_wld, muts_delsel_next, dims=3)
            muts_bensel_wld = cat(muts_bensel_wld, muts_bensel_next, dims=3)
        end
        if occursin("N", data_to_generate)
            muts_delneu_wld = cat(muts_delneu_wld, muts_delneu_next, dims=3)
            muts_benneu_wld = cat(muts_benneu_wld, muts_benneu_next, dims=3)
        end
    end

    wld_stats["x_max_burnin"] = x_max_burnin
    wld_stats["y_max_burnin"] = DEF_Y_MAX
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    wld_stats["n_gens"] = n_gens_total
    
    return Dict("stats"=>wld_stats,"meanf"=>meanf_wld,"pops"=>pops_wld,"delsel"=>muts_delsel_wld,"delneu"=>muts_delneu_wld,
        "bensel"=>muts_bensel_wld,"benneu"=>muts_benneu_wld)
end

# Plotting functions
# ------------------------------------------------

"""
Shows an animated heatmap of ```obj``` from ```gen_start``` to ```gen_end```.

---

```obj```: 3-dimensional array of 2d by-deme data + time 

```gen_start```: start generation

```gen_end```: end generation

```slow_factor```: number of animation frames per generation

```clim```: color bounds (Plots.jl's clim parameter)

```log_base```: if not **-1**, color shows log values with this as base

---

"""
function re_heatmap(obj,gen_start=1,gen_end=DEF_N_GENS_BURNIN+DEF_N_GENS_EXP;slow_factor=1,clim=:default,log_base=-1)
    @gif for i=gen_start:(gen_end*slow_factor-1)
        gen_no = trunc(Int,i/slow_factor)+1
        if all(isnan,obj[:,:,gen_no])
            println("No values found in any deme.")
        end
        if log_base>0 && log_base==1
            heatmap(log.(log_base,obj[:,:,gen_no]'),ylabel="Generation $gen_no",size=(1000,250),clim=clim)
        else
            heatmap(obj[:,:,gen_no]',ylabel="Generation $gen_no",size=(1000,250),clim=clim)
        end
    end
end

"""
Shows population data of ```re``` from ```gen_start``` to ```gen_end```.

---

```re```: range expansion results dictionary

```gen_start```: start generation

```gen_end```: end generation

```slow_factor```: number of animation frames per generation

```clim```: color bounds (Plots.jl's clim parameter)

```log_base```: if not **-1**, color shows log values with this as base

---

"""
function re_heatmap_pops(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["k_capacity"]),log_base=-1)
    re_heatmap(re["pops"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_meanf(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,1),log_base=-1)
    re_heatmap(re["meanf"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end
function re_heatmap_meanf(obj::Array,gen_start,gen_end;slow_factor=1,clim=(0,1),log_base=-1)
    re_heatmap(obj,gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_AAsel(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["AAsel"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_Aasel(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["Aasel"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_aasel(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["aasel"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_AAneu(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["AAneu"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_Aaneu(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["Aaneu"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_aaneu(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["aaneu"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_delsel(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1)
    re_heatmap(re["delsel"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_delneu(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1)
    re_heatmap(re["delneu"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_bensel(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1)
    re_heatmap(re["bensel"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

function re_heatmap_benneu(re,gen_start=1,gen_end=re["stats"]["n_gens"];slow_factor=1,clim=(0,50),log_base=-1)
    re_heatmap(re["benneu"],gen_start,gen_end;slow_factor,clim=clim,log_base=log_base)
end

# Functions pertaining to averaging and expansion front
# ------------------------------------------------

function average_all(data::Array,n_gens::Int;greaterzero=false,divide=true)
    res = Array{Float32}(undef,0)
    for j in 1:n_gens
        push!(res,mean(skipmissing(data[:,:,j])))
    end
    return res
end
function average_all(re::Dict,dataname::String;greaterzero=false,divide=true)
    res = Array{Float32}(undef,0)
    for j in 1:re["stats"]["n_gens"]
        push!(res,mean(filter(!isnan,re[dataname][:,:,j])))
    end
    return res
end

"""
Finds the average value of ```obj``` between all demes at the expansion front of ```re```.

---

```re```: range expansion results dictionary

```obj```: 3-dimensional array of 2d by-deme data + time 

```greaterzero```: if **true**, **>0** values are considered when determining the front (**>=0** values if **false**)

```oneside```: if **true**, approach only from one side (i.e. from the positive direction in axial expansions)

```divide```: if **true**, find average

---

"""
function average_front(re,obj="meanf";greaterzero=false,oneside=false,divide=true)
    average_front(re[obj],re["stats"]["n_gens"],re["stats"]["x_max"],re["stats"]["y_max"];greaterzero=greaterzero,oneside=oneside,divide=divide)
end

function average_front(data_array,n_gens,x_max,y_max;greaterzero=false,oneside=false,divide=true)
    front_array = Array{Float32}(undef,0)
    for j in 1:n_gens
        sum_total = 0
        cnt = 0
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && (isnan(data_array[frontier_x,_y,j]) || (greaterzero && data_array[frontier_x,_y,j] == 0))
                frontier_x -= 1
            end
            if data_array[frontier_x,_y,j]>=0 || (greaterzero && data_array[frontier_x,_y,j]>0)
                sum_total += data_array[frontier_x,_y,j]
                cnt += 1
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && (isnan(data_array[frontier_x,_y,j]) || (greaterzero && data_array[frontier_x,_y,j] == 0))
                    frontier_x += 1
                end
                if data_array[frontier_x,_y,j]>=0 || (greaterzero && data_array[frontier_x,_y,j]>0)
                    sum_total += data_array[frontier_x,_y,j]
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
                while frontier_y != 1 && (isnan(data_array[_x,frontier_y,j]) || (greaterzero && data_array[_x,frontier_y,j] == 0))
                    frontier_y -= 1
                end
                if data_array[_x,frontier_y,j]>=0 || (greaterzero && data_array[_x,frontier_y,j]>0)
                    sum_total += data_array[_x,frontier_y,j]
                    cnt += 1
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && (isnan(data_array[_x,frontier_y,j]) || (greaterzero && data_array[_x,frontier_y,j] == 0))
                    frontier_y += 1
                end
                if data_array[_x,frontier_y,j]>=0 || (greaterzero && data_array[_x,frontier_y,j]>0)
                    sum_total += data_array[_x,frontier_y,j]
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
Finds the front array of ```obj``` in ```re```.

---

```re```: range expansion results dictionary

```obj```: 3-dimensional array of 2d by-deme data + time 

```oneside```: if **true**, approach only from one side (i.e. from the positive direction in axial expansions)

---

"""
function front_array(re,obj="meanf";oneside=false)
    front_array(re[obj],re["stats"]["n_gens"],re["stats"]["x_max"],re["stats"]["y_max"];oneside=oneside)
end

function front_array(data_array,n_gens,x_max,y_max;oneside=false)
    front_arr = fill(NaN,x_max,y_max,n_gens)
    for j in 1:n_gens
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && isnan(data_array[frontier_x,_y,j])
                frontier_x -= 1
            end
            if !isnan(data_array[frontier_x,_y,j])
                front_arr[frontier_x,_y,j]=data_array[frontier_x,_y,j]
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && isnan(data_array[frontier_x,_y,j])
                    frontier_x += 1
                end
                if !isnan(data_array[frontier_x,_y,j])
                    front_arr[frontier_x,_y,j]=data_array[frontier_x,_y,j]
                end
            end
        end

        if !oneside
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && isnan(data_array[_x,frontier_y,j])
                    frontier_y -= 1
                end
                if !isnan(data_array[_x,frontier_y,j])
                    front_arr[_x,frontier_y,j]=data_array[_x,frontier_y,j]
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && isnan(data_array[_x,frontier_y,j])
                    frontier_y += 1
                end
                if !isnan(data_array[_x,frontier_y,j]>0)
                    front_arr[_x,frontier_y,j]=data_array[_x,frontier_y,j]
                end
            end
        end
    end
    return front_arr
end

# mean front fitness (or other data)
"""
Normalises ```obj``` in ```re``` using the "maximum normalisation" method: after the last burn-in generation, divide by the maximum of each generation.

---

```re```: range expansion results dictionary

```obj```: 3-dimensional array of 2d by-deme data + time 

---

"""
function norm_maximum(re,obj="meanf")
    norm_maximum(re[obj],re["stats"]["n_gens_burnin"],re["stats"]["n_gens_exp"])
end

function norm_maximum(data_array,n_gens_burnin,n_gens_exp)
    normal_array = copy(data_array)
    for j in 1:n_gens_exp
        gen_max = maximum(data_array[:,:,n_gens_burnin+j])
        normal_array[:,:,n_gens_burnin+j] /= gen_max
    end
    return normal_array
end

"""
Normalises ```obj``` in ```re``` using the "maximum normalisation" method: after the last burn-in generation, divide by the constant value of average fitness over all demes at the (last burn-in generation+1)=onset generation

---

```re```: range expansion results dictionary

```obj```: 3-dimensional array of 2d by-deme data + time 

```offset``` - offset from the onset generation

---

"""
function norm_onset_mean(re::Dict,obj::String="meanf",offset=0)
    norm_onset_mean(re[obj],re["stats"]["n_gens_burnin"],offset=offset)
end

function norm_onset_mean(data_array::Array,n_gens_burnin::Int,offset=0)
    normal_array = copy(data_array)

    sum = 0
    count = 0
    for u in data_array[:,:,n_gens_burnin+1+offset]
        if u > 0
            sum += u
            count += 1
        end
    end
    gen_average = sum/count

    normal_array[:,:,n_gens_burnin+1:end] /= gen_average
    return normal_array
end

#= function normalise_front_by_onset_mean(average_1d_array)
    normal_array = copy(average_1d_array)
    normal_array[BURN_IN_GEN_N+1:end] /= average_1d_array[BURN_IN_GEN_N+1]
    return normal_array
end

function normalise_front_by_max(average_1d_array,meanf_array)
    normal_array = copy(average_1d_array)
    for j in 1:(TOTAL_GEN_N-BURN_IN_GEN_N)
        gen_max = maximum(meanf_array[:,:,BURN_IN_GEN_N+j])
        normal_array[BURN_IN_GEN_N+j] /= gen_max
    end
    return normal_array
end

function find_front_array_muts(data_array,muts_array;oneside=false)
    res_muts = zeros(Float32,y_max,n_gen)
    for j in 1:n_gen
        # scanning every y: side 1

        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && data_array[frontier_x,_y,j] < 0
                frontier_x -= 1
            end
            if data_array[frontier_x,_y,j]>0
                res_muts[_y,j]=muts_array[frontier_x,_y,j]
            end
        end
        # scanning every y: side 2
        # add later
    end
    return res_muts
end =#

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

println("2d.jl successfully loaded.")