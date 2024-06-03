using StatsBase, Distributions, Random, Plots, SpecialFunctions, Serialization

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

# Gene parameters
DEF_N_LOCI = 1000
DEF_N_SEL_LOCI = 500
DEF_MUT_RATE = 0.7567 # genome-wide
DEF_MIGR_RATE = 0.1
DEF_S_SEL_COEF = 0.002
DEF_H_DOMIN_COEF = 0
DEF_PROP_OF_DEL_MUTS = 0.9
MIGR_DIRS_8 = [[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]
MIGR_DIRS_6 = [[-1,0],[0,-1],[-1,1],[0,1],[1,0],[1,1]]
MIGR_DIRS_4 = [[-1,0],[0,-1],[0,1],[1,0]]
DEF_MIGR_DIRS = MIGR_DIRS_4

# Expansion parameters
DEF_X_MAX_BURNIN = 5
DEF_N_GENS_BURNIN = 10
DEF_X_MAX_EXP = DEF_X_MAX
DEF_Y_MAX_EXP = DEF_Y_MAX
DEF_N_GENS_EXP = 20
DEF_MIGR_MODE = 0
DEF_DATA_TO_GENERATE = "FP"