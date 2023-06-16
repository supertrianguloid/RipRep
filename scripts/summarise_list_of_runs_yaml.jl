#!julia
using YAML
include(Base.source_dir()*"/../src/parser.jl")
include(Base.source_dir()*"/../src/spectroscopy.jl")
include(Base.source_dir()*"/../src/utilities.jl")

OUTPUT_DIRECTORY = "/home/lbowes/ANALYSIS/"
THERM = 1000
BINSIZE = 10

ensembles = Dict()

for line in readlines(ARGS[1])
    ensemble_path = OUTPUT_DIRECTORY * replace(line, "/" => "_")[2:end] * "/"
    ensure_directory_exists(ensemble_path)
    println("Processing " * line)
    ens = load_ensemble(line)
    ensembles[line] = ens.global_metadata
    if ens.global_metadata[:nconfs] > 1100
	thermalise!(ens, THERM)
	ensembles[line][:PCAC] = bootstrap_effective_pcac(ens.analysis, BINSIZE)
    end
    YAML.write_file(ensemble_path * "analysis.yml", ensembles[line])
end

