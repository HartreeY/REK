# fitness-multiple-dims
A set of programs in Julia designed to efficiently simulate range expansions and study their genetics and population dynamics.
This set of programs has been used in the study "The evolution of fitness during range expansions in multiple dimensions". You can find the preprint at [[https://www.biorxiv.org/content/10.1101/2023.12.29.573608v2]].

# Main use cases
To begin using this set of tools, start by install the required package, which should take only around a minute. Access the *programs* folder and run the script *init.jl*.

After installing the required packages, select an appropriate script according to the dimensionality of your required simulation. For example, the *2d.ipynb* script includes a handful of commands for simulating range expansions under different conditions. To run a simple simulation, try running the *rangeexp_axial()* function:
```
test, test_stats = rangeexp_axial()
```
