# Range Expansions Simulation Kit (RESK)
A set of programs in Julia designed to efficiently simulate range expansions and study their genetics and population dynamics.
This set of programs has been used in the study "The evolution of fitness during range expansions in multiple dimensions". You can find the preprint at [[https://www.biorxiv.org/content/10.1101/2023.12.29.573608v2]].

# How to use
## Prerequisites
To begin using this set of tools, you need to have Julia 1.9+ installed, along with the required packages. To install the packages, please run the initialisation script: access the *programs* folder and run the script *init.jl*. This should only take around 2 minutes.

Once you have the required packages, you can use the annotated methods of this package on your own, or you can follow one of several examples in the *programs* folder.

## Main use cases
Include the appropriate script according to the dimensionality of your required simulation and use its methods. For example, the *2d.jl* script includes a handful of commands for simulating range expansions in two dimensions.

To run a simulation once, use one of the following methods: **rangeexp_axial**, **rangeexp_radial**. These output a fixed dictionary that includes statistics and expansion data. The types of data within it are determined by the *data_to_generate* argument. It can take on the following values:
- **F** - **meanf** (deme-average fitness)
- **P** - **pops** (deme populations)
- **S** - **AAsel**, **Aasel** and **aasel** (deme-average number of homo- [**AA**, **aa**] and heterozygous [**Aa**] selected loci)
- **N** - **AAneu**, **Aaneu** and **aaneu** (deme-average number of homo- and heterozygous neutral loci)
The above can be combined and should be passed in a string. For example,
```
test = rangeexp_axial(15, 30; data_to_generate="SF", y_max=5)
```
will be a
```
Dict{String, Any} with 9 entries:
  "AAsel" => Float32[0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.…
  "Aasel" => Float32[0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.…
  "aasel" => Float32[0.0 0.0 … 0.0 0.0; 25.0 25.0 … 0.0 25.0; … ; 0.0 0.0 … 0.0…
  "pops"  => NaN
  "meanf" => Float32[-1.0 -1.0 … -1.0 -1.0; 0.05 0.05 … -1.0 0.05; … ; -1.0 -1.…
  "aaneu" => NaN
  "Aaneu" => NaN
  "AAneu" => NaN
  "stats" => Dict{String, Any}("y_max_burnin"=>10, "x_max"=>100, "migr_mode"=>[…
```

This expansion data can be plotted and worked with. To plot expansion data, use the unique functions that start with *re_*. For example,
```
re_heatmap_AAsel(test; log_factor=1.02)
```
will output the average number of selected homozygous mutant loci for in a deme:
![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme0.gif?raw=true)

Here's an example of a longer axial simulation in 2D:
```
test = rangeexp_axial_inf(100,1000;data_to_generate="FPSN",prop_of_sel_loci=0.8,y_max=8,migr_mode="diag1/2")
```
![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme1.gif?raw=true)