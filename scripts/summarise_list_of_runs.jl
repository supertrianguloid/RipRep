#!julia
include(Base.source_dir()*"/../src/parser.jl")

ensembles = []

for line in readlines(ARGS[1])
    println("Processing " * line)
    push!(ensembles, load_ensemble(line).global_metadata)
end

println("|beta|csw|m0|T|L|nconfs|")
for ensemble in ensembles
    println("|" * string(ensemble[:β]) * "|" * string(ensemble[:csw]) * "|" * string(ensemble[:m0]) * "|" * string(ensemble[:geometry][0]) * "|" * string(ensemble[:geometry][1]) * "|" * string(ensemble[:nconfs]) * "|")
end
