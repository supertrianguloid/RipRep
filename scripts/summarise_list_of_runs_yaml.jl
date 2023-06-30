#!julia
using YAML
using Dates
using LoggingExtras

global_logger(MinLevelLogger(FileLogger(ARGS[2]), Logging.Info))



include(Base.source_dir()*"/../src/parser.jl")
include(Base.source_dir()*"/../src/spectroscopy.jl")
include(Base.source_dir()*"/../src/wilson.jl")
include(Base.source_dir()*"/../src/utilities.jl")

#Turn off plot display because we are headless
ENV["GKSwstype"]="nul"

OUTPUT_DIRECTORY = "/home/lbowes/ANALYSIS/" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS") * "/"
THERM = 1000
BINSIZE = 10
NO_FIT_POINTS = 4

WF_REF = 1.0

ensembles = Dict()

ensure_directory_exists(OUTPUT_DIRECTORY)

list_of_ensembles = YAML.load_file(ARGS[1])

function process_ensemble(line, ensemble_data)
    ensemble_path = OUTPUT_DIRECTORY * replace(line, "/" => "_")[2:end] * "/"
    ensure_directory_exists(ensemble_path)
    with_logger(MinLevelLogger(FileLogger(ensemble_path * "analysis.log"), Logging.Info)) do
        wf = nothing
        try
            wf = ensemble_data["wf"]
            wf = load_wilsonflow(wf)
            @info "Successfully loaded Wilson flow data"
        catch e
        end

        save_figure(name) = savefig(ensemble_path * name)
        ens = load_ensemble(line)
        analysis = deepcopy(ens.global_metadata)
        contains_nans = ens.global_metadata[:nan_confs] != []
        if ens.global_metadata[:nconfs] > 1100 && !contains_nans
            corrs = [:g5_folded, :gk_folded, :id_folded]
            T = ens.global_metadata[:geometry][0]
            L = ens.global_metadata[:geometry][1]
            T_middle = T ÷ 2
            fit_window = (T_middle - NO_FIT_POINTS):T_middle
            thermalise!(ens, THERM)
            try
                @info "Plotting plaquette..."
                plot_plaquette(ens)
                save_figure("plaquette.pdf")
            catch e
                @error "Failed!"
            end
            try
                @info "PCAC mass..."
                analysis[:pcac] = bootstrap_effective_pcac(ens.analysis, BINSIZE)
                plot_pcac_mass(ens, BINSIZE)
                save_figure("pcac_mass.pdf")
                analysis[:m_pcac] = fit_pcac_mass(ens, BINSIZE, fit_window)
                plot_pcac_fit(ens, BINSIZE, fit_window, 1:T_middle)
                save_figure("pcac_mass_fit.pdf")
            catch e
                @error "Failed!"
            end
            try
                @info "Fps..."
                analysis[:effective_fps] = bootstrap_effective_fps(ens.analysis, T, L, BINSIZE)
                plot_fps(ens, BINSIZE)
                save_figure("fps.pdf")
                analysis[:fps] = fit_fps(ens, BINSIZE, fit_window)
                plot_fps_fit(ens, BINSIZE, fit_window, 1:T_middle)
                save_figure("fps_fit.pdf")
            catch e
                @error "Failed!"
            end
            try
                @info "Gps..."
                analysis[:effective_gps] = bootstrap_effective_gps(ens.analysis, ens.global_metadata[:geometry][1], BINSIZE)
                plot_gps(ens, BINSIZE)
                save_figure("gps.pdf")
                analysis[:gps] = fit_gps(ens, BINSIZE, fit_window)
                plot_gps_fit(ens, BINSIZE, fit_window, 1:T_middle)
                save_figure("gps_fit.pdf")
            catch e
                @error "Failed!"
            end
            for corr in corrs
                try
                    @info String(corr) * "..."
                    analysis[Symbol("effective_"* String(corr))] = bootstrap_effective_mass(ens.analysis, corr, BINSIZE)
                    plot_effective_mass(ens, corr, BINSIZE)
                    save_figure("effective_mass_" * String(corr) * ".pdf")
                    analysis[corr] = fit_effective_mass(ens, corr, BINSIZE, fit_window)
                    plot_effective_mass_fit(ens, corr, BINSIZE, fit_window, 1:T_middle)
                    save_figure("effective_mass_" * String(corr) * "_fit.pdf")
                catch e
                    @error "Failed!"
                end
            end
            try
                @info "mv/mpi..."
                analysis[:effective_ratio_mv_mpi] = bootstrap_effective_mass_ratio(ens.analysis, :gk_folded, :g5_folded, BINSIZE)
                plot_effective_mass_ratio(ens, :gk_folded, :g5_folded, BINSIZE)
                save_figure("effective_ratio_mv_mpi.pdf")
                analysis[:ratio_mv_mpi] = [last(analysis[:effective_mass_ratio_mv_mpi][1]), last(analysis[:effective_mass_ratio_mv_mpi][2])]
            catch e
                @error "Failed!"
            end
            if wf != nothing
                @info "Wilson flow..."
                try
                    @info "t²E..."
                    plot_t2e(wf, binsize = BINSIZE)
                    save_figure("t2e.pdf")
                catch e
                    @error "Failed!"
                end
                try
                    @info "W..."
                    plot_w(wf, binsize = BINSIZE)
                    save_figure("W.pdf")
                catch e
                    @error "Failed!"
                end
                try
                    @info "Topological Charge..."
                    plot_tc(wf)
                    save_figure("tc.pdf")
                catch e
                    @error "Failed!"
                end
                try
                    @info "Topological Charge Histogram..."
                    plot_tc_hist(wf)
                    save_figure("tc_hist.pdf")
                catch e
                    @error "Failed!"
                end
                try
                    @info "Calculating w0..."
                    w0 = auto_w0(wf)
                    analysis[:w0] = w0
                catch e
                    @error "Failed!"
                end
                try
                    @info "Calculating mpiw0..."
                    analysis[:mpiw0] = propagate_product(analysis[:w0], analysis[:g5_folded])
                catch e
                    @error "Failed!"
                end
            end
            try
                @info "Calculating mpi^2..."
                analysis[:mpi2] = propagate_product(analysis[:g5_folded], analysis[:g5_folded])
            catch e
                @error "Failed!"
            end
            try
                @info "Calculating mpiL..."
                analysis[:mpiL] = analysis[:g5_folded][1]*L
            catch e
                @error "Failed!"
            end
        end
        
    end
    YAML.write_file(ensemble_path * "analysis.yml", analysis)
    if contains_nans
        touch(ensemble_path * "CONTAINS_NANS")
    end
    
end

for line in [k for k in keys(list_of_ensembles)]
    @info "Processing " * line
    analysis = process_ensemble(line, list_of_ensembles[line])
    ensemble[line] = analysis
end


YAML.write_file(OUTPUT_DIRECTORY * "ensembles.yml", ensembles)
