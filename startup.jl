# automatic load of things related to all projects go here

import Pkg  # or using Pkg

current_directory =  @__DIR__() 
print( "Current directory is: ", current_directory, "\n\n" )
Pkg.activate(source_directory)  # so now you activate the package
Base.active_project()  
push!(LOAD_PATH, source_directory)  # add the directory to the load path, so it can be found

pkgs_startup = [  
    "Revise", "Logging", "OhMyREPL",
    "Setfield", "Memoization",
    "DataFrames", "CSV", "RData",    
    "StatsBase", "Statistics", "Distributions", "Random", 
    "PlotThemes", "Colors", "ColorSchemes", 
    "Plots", "StatsPlots",
    "MKL", "ForwardDiff", "PDMats",   
    "MultivariateStats",  "StatsModels", "StatsFuns",
    "StaticArrays", "LazyArrays", "FillArrays", "SparseArrays", "Graphs",
    "Flux", "Optim", "KernelFunctions", "AbstractGPs", 
    "Distributions", "DistributionsAD", "AdvancedVI",
    "Interpolations", "LinearAlgebra",
    "Turing"  # should be last
]
  
if @isdefined pkgs 
    pkgs = unique!( [pkgs_startup; pkgs] )
else 
    pkgs = pkgs_startup
end
    

print( "Loading libraries:\n\n" )
for pk in pkgs; 
    print(pk, "\n"); @eval using $(Symbol(pk)); 
end    
 

colorscheme!("Monokai24bit") # for REPL

# some convenience functions
print( "\n\n", "Additional functions: \n\n" )

showall(x) = show(stdout, "text/plain", x)
Base.dim = Base.size  # 
 