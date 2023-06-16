#!julia
using YAML
include(Base.source_dir()*"/../src/parser.jl")
include(Base.source_dir()*"/../src/spectroscopy.jl")
include(Base.source_dir()*"/../src/utilities.jl")

#Turn off plot display because we are headless
ENV["GKSwstype"]="nul"

OUTPUT_DIRECTORY = "/home/lbowes/ANALYSIS/"
THERM = 1000
BINSIZE = 10

ensembles = Dict()

for line in readlines(ARGS[1])
    ensemble_path = OUTPUT_DIRECTORY * replace(line, "/" => "_")[2:end] * "/"
    ensure_directory_exists(ensemble_path)
    println("Processing " * line)
    save_figure(name) = savefig(ensemble_path * name)
    ens = load_ensemble(line)
    ensembles[line] = ens.global_metadata
    contains_nans = ens.global_metadata[:nan_confs] != []
    if ens.global_metadata[:nconfs] > 1100 && !contains_nans
        corrs = [:g5_folded, :gk_folded, :id_folded]
	thermalise!(ens, THERM)
        try
            plot_plaquette(ens)
            save_figure("plaquette.pdf")
        catch e
            @error "Failed!"
        end
        try
            @info "PCAC mass..."
            ensembles[line][:pcac] = bootstrap_effective_pcac(ens.analysis, BINSIZE)
            plot_pcac_mass(ens, BINSIZE)
            save_figure("pcac_mass.pdf")
        catch e
            @error "Failed!"
        end
        try
            @info "Fps..."
            ensembles[line][:fps] = bootstrap_effective_fps(ens.analysis, ens.global_metadata[:geometry][0], ens.global_metadata[:geometry][1], BINSIZE)
            plot_fps(ens, BINSIZE)
            save_figure("fps.pdf")
        catch e
            @error "Failed!"
        end
        try
            @info "Gps..."
            ensembles[line][:gps] = bootstrap_effective_gps(ens.analysis, ens.global_metadata[:geometry][1], BINSIZE)
            plot_gps(ens, BINSIZE)
            save_figure("gps.pdf")
        catch e
            @error "Failed!"
        end
        for corr in corrs
            try
                @info String(corr) * "..."
                ensembles[line][corr] = bootstrap_effective_mass(ens.analysis, corr, BINSIZE)
                plot_effective_mass(ens, corr, BINSIZE)
                save_figure("effective_mass_" * String(corr) * ".pdf")
            catch e
                @error "Failed!"
            end
        end
    end
    YAML.write_file(ensemble_path * "analysis.yml", ensembles[line])
    if contains_nans
        touch(ensemble_path * "CONTAINS_NANS")
    end
end

