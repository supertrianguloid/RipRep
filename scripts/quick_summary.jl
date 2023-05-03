#!/usr/bin/julia
include(Base.source_dir()*"/../src/parser.jl")
print(load_ensemble(ARGS[1]).global_metadata)
