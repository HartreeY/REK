include("core.jl")

# Simulation functions
# ------------------------------------------------

"""
Calculates average fitness and average hetero- and homozygosities in a deme.
"""
function muts_by_sel_neu(inds,s_sel_coef,n_loci,sel_loci=[])
    len = length(inds)
    muts_sel = 0
    muts_neu = 0
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
        muts1s += muts_AA_sel
        muts2s += muts_Aa_sel
        muts1ns += muts_AA_neu
        muts2ns += muts_Aa_neu
        muts3s += length(sel_loci) - muts_AA_sel - muts_Aa_sel
        muts3ns += n_loci - length(sel_loci) - muts_AA_neu - muts_Aa_neu
    end

    muts_sel /= len
    muts_neu /= len
    return muts_sel,muts_neu,fits
end

@inbounds function mutate(person,mut_rate,n_loci,s_sel_coef)
    get_mutation_random = rand(Poisson(mut_rate))
    @fastmath @inbounds for _ in 1:get_mutation_random
        cnt_del_muts = 0
        cnt_ben_muts = 0
        pos_alter = sample(1:n_loci*2)

        if rand() < mut_rate
            person[pos_alter] *= 1 - s_sel_coef
            cnt_del_muts += 1
        else
            person[pos_alter] *= 1 + s_sel_coef
            cnt_ben_muts += 1
        end

        return cnt_del_muts, cnt_ben_muts
    end
end

@inbounds function crossover(person,n_loci)
    for i in 1:n_loci
        lr = rand(1:2)
        person[i] = lr==1 ? person[i] : person[i+n_loci]
    end
end

@inbounds function mate(person1,person2,n_loci)
    lr1 = rand(1:2)==1 ? (1:n_loci) : ((n_loci+1):(n_loci*2))
    lr2 = rand(1:2)==1 ? (1:n_loci) : ((n_loci+1):(n_loci*2))
    return vcat(person1[lr1],person2[lr2])
end

@inbounds function build_next_gen(pnt_wld,pnt_wld_stats,pnt_meanf_wld=NaN,pnt_pops_wld=NaN,pnt_muts_wld=NaN;
    x_max_migr=NaN,y_max_migr=NaN,migr_mode=0,x_bottleneck=NaN,refl_walls=false)

    # Determine the number of offspring for each deme
    next_gen_posits = []
    next_gen_pops = fill(NaN,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
    for x in 1:pnt_wld_stats["x_max"],y in 1:pnt_wld_stats["y_max"]
        if isassigned(pnt_wld,x,y) && length(pnt_wld[x,y])>0
            n_ppl_at_deme = length(pnt_wld[x,y])
            expected_offspring = n_ppl_at_deme * (pnt_wld_stats["r_prolif_rate"]/(1 + (n_ppl_at_deme*(pnt_wld_stats["r_prolif_rate"]-1))/pnt_wld_stats["k_capacity"]))
            next_gen_pops[x,y] =  rand(Poisson(expected_offspring))
            if next_gen_pops[x,y]>0
                push!(next_gen_posits,[x,y])
            end
        end
    end
    
    # Define the habitat (world) and the data arrays in the next generation
    wld_next = Array{Array{Array{Float32}}}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
    mean_fitn_next = NaN
    pops_next = NaN
    muts_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    meanf_out = false
    pops_out = false
    muts_out = false
    if pnt_meanf_wld isa Array{Float32, 3}
        meanf_out = true
        mean_fitn_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(mean_fitn_next,NaN)
    end
    if pnt_pops_wld isa Array{Int32, 3}
        pops_out = true
        pops_next = zeros(Int32,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
    end
    if pnt_muts_wld isa Array{Float32, 3}
        muts_out = true
        muts_next = Array{Float32}(undef,pnt_wld_stats["x_max"],pnt_wld_stats["y_max"])
        fill!(muts_next,NaN)
    end
    

    for deme in next_gen_posits
        inds_at_pos = pnt_wld[deme...]

        fitns = []
        cnt_muts,fitns =
            muts_by_sel_neu(inds_at_pos,pnt_wld_stats["s_sel_coef"],pnt_wld_stats["n_loci"],pnt_wld_stats["sel_loci"])
        
        sum_fitn = sum(fitns)
        fitns /= sum_fitn

        if meanf_out
            mean_fitn_next[deme...] = mean(fitns)
        end
        if muts_out
            muts_next[deme...] = cnt_muts
        end

        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            birth_count = 0
            for _ in 1:next_generation_size
                mom = wsample(inds_at_pos,fitns)
                dad = wsample(inds_at_pos,fitns)

                gamete_mom = copy(mom)
                gamete_dad = copy(dad)

                crossover(gamete_mom,pnt_wld_stats["n_loci"])
                crossover(gamete_dad,pnt_wld_stats["n_loci"])
                mutate(gamete_mom,pnt_wld_stats["mut_rate"],pnt_wld_stats["n_loci"],pnt_wld_stats["s_sel_coef"])
                mutate(gamete_dad,pnt_wld_stats["mut_rate"],pnt_wld_stats["n_loci"],pnt_wld_stats["s_sel_coef"])
                mate_result = mate(gamete_mom,gamete_dad,pnt_wld_stats["n_loci"])

                p_lat,p_diag = get_migr_params(migr_mode)

                move_x = 0
                move_y = 0
                migr_res = rand()
                if rand()<pnt_wld_stats["migr_rate"] && migr_res < p_lat+p_diag
                    if migr_res < p_lat
                        dir = sample(MIGR_DIRS_4)
                    elseif migr_res < p_lat+p_diag
                        dir = sample(MIGR_DIRS_8)
                    end

                    # Raw migration results
                    move_x = dir[1]
                    move_y = dir[2]

                    # Nullify migration on certain conditions
                    if !isnan(x_bottleneck) && (deme[1]+move_x==x_bottleneck && deme[2]+move_y!=ceil(pnt_wld_stats["y_max"]/2)) # bottleneck barrier check
                        move_x = 0
                        move_y = 0
                    else
                        if deme[1]+move_x > x_max_migr || deme[1]+move_x < 1 # burn-in area check
                            move_x = refl_walls ? -move_x : 0
                        end
                        if deme[2]+move_y > pnt_wld_stats["y_max"] || deme[2]+move_y < 1
                            move_y = refl_walls ? -move_y : 0
                        end
                    end
                end

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
    
    return wld_next, mean_fitn_next, pops_next, muts_next
end

function create_empty_world(x_max=DEF_X_MAX,y_max=DEF_Y_MAX;name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"),k_capacity=DEF_K_CAPACITY,
    r_prolif_rate=DEF_R_PROLIF_RATE,n_loci=DEF_N_LOCI,n_sel_loci=DEF_N_SEL_LOCI,
    mut_rate=DEF_MUT_RATE,migr_rate=DEF_MIGR_RATE,migr_dirs=DEF_MIGR_DIRS,s_sel_coef=DEF_S_SEL_COEF,prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

    wld = Array{Array{Array{Float32}}}(undef,x_max,y_max) # array of fitness values of all individuals in space

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
        "migr_dirs" => migr_dirs,
        "s_sel_coef" => s_sel_coef,
        "prop_of_del_muts" => prop_of_del_muts

        #"rangeexps" => []
    )

    return wld,wld_stats
end

function fill_random_demes(pnt_wld,pnt_wld_stats,x_max_fill,y_max_fill,n_demes_to_fill=DEF_N_DEMES_STARTFILL)

    possible_init_coords = [collect(x) for x in Iterators.product(1:x_max_fill, 1:y_max_fill)]
    init_coords = sample(possible_init_coords,n_demes_to_fill;replace=false)
    pnt_wld_stats["sel_loci"] = randperm(pnt_wld_stats["n_loci"])[1:pnt_wld_stats["n_sel_loci"]]

    for coord in init_coords
        if !isassigned(pnt_wld,coord...)
            pnt_wld[coord...] = []
        end
        for _ in 1:pnt_wld_stats["k_capacity"]
            push!(pnt_wld[coord...],ones(pnt_wld_stats["n_loci"]*2))
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

function rangeexp_axial(n_gens_burnin=DEF_N_GENS_BURNIN,n_gens_exp=DEF_N_GENS_EXP;x_max_burnin=DEF_X_MAX_BURNIN,x_max_exp=DEF_X_MAX_EXP,y_max=DEF_Y_MAX,migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE,wld=NaN,wld_stats=NaN,name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS")) # expansion along x model 1

    meanf_wld = NaN
    pops_wld = NaN
    muts_sel_wld = NaN
    muts_neu_wld = NaN

    if !(wld isa Array{Float32, 3})
        #println("No world provided. Creating a new world.")
        wld,wld_stats = create_empty_world(x_max_exp,y_max;name=name)
        fill_random_demes(wld,wld_stats,Int(x_max_exp/20),y_max)
    end

    if occursin("F", data_to_generate)
        meanf_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("P", data_to_generate)
        pops_wld = Array{Int32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("S", data_to_generate)
        muts_AAsel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_Aasel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
        muts_aasel_wld = Array{Float32}(undef, wld_stats["x_max"], wld_stats["y_max"], 0)
    end
    if occursin("M", data_to_generate)
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
        if occursin("M", data_to_generate)
            muts_AAneu_wld = cat(muts_AAneu_wld, muts_AAneu_next, dims=3)
            muts_Aaneu_wld = cat(muts_Aaneu_wld, muts_Aaneu_next, dims=3)
            muts_aaneu_wld = cat(muts_aaneu_wld, muts_aaneu_next, dims=3)
        end
    end

#=  append!(wld_stats["rangeexps"],Dict(
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
    
#=     res = []
    function needed_data(symb,arr,res)
        if occursin("P", data_to_generate)
            push!(res,Ref(arr))
        end
    end
    needed_data("P",pops_wld) =#
    return Dict("stats"=>wld_stats,"meanf"=>meanf_wld,"pops"=>pops_wld,"AAsel"=>muts_AAsel_wld,"Aasel"=>muts_Aasel_wld,
        "aasel"=>muts_aasel_wld,"AAneu"=>muts_AAneu_wld,"Aaneu"=>muts_Aaneu_wld,"aaneu"=>muts_aaneu_wld)
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
function re_heatmap(obj,gen_start=1,gen_end=DEF_N_GENS_BURNIN+DEF_N_GENS_EXP,slow_factor=1;clim=:default,log_base=-1)
    @gif for i=gen_start:(gen_end*slow_factor-1)
        gen_no = trunc(Int,i/slow_factor)+1
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
function re_heatmap_pops(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,re["stats"]["k_capacity"]),log_base=-1)
    re_heatmap(re["pops"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_meanf(re::Dict,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,1),log_base=-1)
    re_heatmap(re["meanf"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end
function re_heatmap_meanf(obj::Array,gen_start,gen_end,slow_factor=1;clim=(0,1),log_base=-1)
    re_heatmap(obj,gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_AAsel(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["AAsel"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_Aasel(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["Aasel"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_aasel(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["aasel"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_AAneu(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["AAneu"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_Aaneu(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["Aaneu"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

function re_heatmap_aaneu(re,gen_start=1,gen_end=re["stats"]["n_gens"],slow_factor=1;clim=(0,re["stats"]["n_loci"]-length(re["stats"]["sel_loci"])),log_base=-1)
    re_heatmap(re["aaneu"],gen_start,gen_end,slow_factor;clim=clim,log_base=log_base)
end

# Functions pertaining to averaging and expansion front
# ------------------------------------------------

function average_all(data::Array,n_gens::Int;leqzero=false,divide=true)
    res = Array{Float64}(undef,0)
    for j in 1:n_gens
        push!(res,mean(skipmissing(data[:,:,j])))
    end
    return res
end
function average_all(re::Dict,dataname::String;leqzero=false,divide=true)
    res = Array{Float64}(undef,0)
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

```leqzero```: if **true**, approach only from one side (i.e. from the positive direction in axial expansions)

```oneside```: if **true**, approach only from one side (i.e. from the positive direction in axial expansions)

```divide```: if **true**, find average

---

"""
function average_front(re,obj="meanf";leqzero=false,oneside=false,divide=true)
    average_front(re[obj],re["stats"]["n_gens"],re["stats"]["x_max"],re["stats"]["y_max"];leqzero=leqzero,oneside=oneside,divide=divide)
end

function average_front(data_array,n_gens,x_max,y_max;leqzero=false,oneside=false,divide=true)
    front_array = Array{Float64}(undef,0)
    for j in 1:n_gens
        sum_total = 0
        cnt = 0
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && ((leqzero && data_array[frontier_x,_y,j] <= 0) || (!leqzero && data_array[frontier_x,_y,j] < 0))
                frontier_x -= 1
            end
            if data_array[frontier_x,_y,j]>0
                sum_total += data_array[frontier_x,_y,j]
                cnt += 1
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && ((leqzero && data_array[frontier_x,_y,j] <= 0) || (!leqzero && data_array[frontier_x,_y,j] < 0))
                    frontier_x += 1
                end
                if data_array[frontier_x,_y,j]>0
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
                while frontier_y != 1 && ((leqzero &&  data_array[_x,frontier_y,j] <= 0) || (!leqzero && data_array[_x,frontier_y,j] < 0))
                    frontier_y -= 1
                end
                if data_array[_x,frontier_y,j]>0
                    sum_total += data_array[_x,frontier_y,j]
                    cnt += 1
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && ((leqzero && data_array[_x,frontier_y,j] <= 0) || (!leqzero && data_array[_x,frontier_y,j] < 0))
                    frontier_y += 1
                end
                if data_array[_x,frontier_y,j]>0
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
    front_arr = zeros(x_max,y_max,n_gens)
    for j in 1:n_gens
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && data_array[frontier_x,_y,j] < 0
                frontier_x -= 1
            end
            if data_array[frontier_x,_y,j]>0
                front_arr[frontier_x,_y,j]=data_array[frontier_x,_y,j]
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && data_array[frontier_x,_y,j] < 0
                    frontier_x += 1
                end
                if data_array[frontier_x,_y,j]>0
                    front_arr[frontier_x,_y,j]=data_array[frontier_x,_y,j]
                end
            end
        end

        if !oneside
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && data_array[_x,frontier_y,j] < 0
                    frontier_y -= 1
                end
                if data_array[_x,frontier_y,j]>0
                    front_arr[_x,frontier_y,j]=data_array[_x,frontier_y,j]
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && data_array[_x,frontier_y,j] < 0
                    frontier_y += 1
                end
                if data_array[_x,frontier_y,j]>0
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
    plot!(x_range*x_scale_factor,t2,label=lbl2)
end

function re_plot_avrelselneu!(re::Dict,dataname::String,x_range=(1:Int(re["stats"]["x_max"]));x_scale_factor=1,sel=true,overlay=false)
    re_plot_avrelselneu(re,dataname,x_range;x_scale_factor=x_scale_factor,sel=sel,overlay=true)
end

println("2d.jl successfully loaded.")