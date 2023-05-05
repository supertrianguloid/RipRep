#!julia
include(Base.source_dir()*"/../src/parser.jl")
using CSV

ensembles = []

for line in readlines(ARGS[1])
    push!(ensembles, only(load_ensemble(line).global_metadata))
end

CSV.write(ARGS[2], DataFrame(ensembles))
