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
    contains_nans = ens.global_metadata[:nan_confs] != []
    if ens.global_metadata[:nconfs] > 1100 && !contains_nans
	thermalise!(ens, THERM)
	ensembles[line][:pcac] = bootstrap_effective_pcac(ens.analysis, BINSIZE)
	ensembles[line][:fps] = bootstrap_effective_fps(ens.analysis, ens.global_metadata[:geometry][0], ens.global_metadata[:geometry][1], BINSIZE)
	ensembles[line][:gps] = bootstrap_effective_gps(ens.analysis, ens.global_metadata[:geometry][1], BINSIZE)
	ensembles[line][:g5] = bootstrap_effective_mass(ens.analysis, :g5_folded, BINSIZE)
	ensembles[line][:gk] = bootstrap_effective_mass(ens.analysis, :gk_folded, BINSIZE)
	ensembles[line][:id] = bootstrap_effective_mass(ens.analysis, :id_folded, BINSIZE)
    end
    YAML.write_file(ensemble_path * "analysis.yml", ensembles[line])
    if contains_nans
        touch(ensemble_path * "CONTAINS_NANS")
    end
end

