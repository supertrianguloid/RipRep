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
DEFAULT_BINSIZE = 10
NO_FIT_POINTS = 4
DEFAULT_TUNE_BINSIZES = collect(1:6:80)

WF_REF = 1.0

ensembles = Dict()

ensure_directory_exists(OUTPUT_DIRECTORY)

list_of_ensembles = YAML.load_file(ARGS[1])

function get_key_or_nothing(dict, key)
    try
	return dict[key]
    catch
        return nothing
    end
end

function get_binsize_tune(dict, key)
    bs = get_key_or_nothing(dict, key)
    if bs == nothing
        return DEFAULT_BINSIZE, true
    end
    return bs, false
end
        

function process_ensemble(line, ensemble_data)
    ensemble_path = OUTPUT_DIRECTORY * replace(line, "/" => "_")[2:end] * "/"
    ensure_directory_exists(ensemble_path)
    with_logger(FileLogger(ensemble_path * "analysis.log")) do
        wf = nothing
        try
            wf = ensemble_data["wf"]
            wf = load_wilsonflow(wf)
            @info "Successfully loaded Wilson flow data"
        catch e
        end
        function save_figure(name)
            savefig(ensemble_path * name * ".pdf")
            savefig(ensemble_path * name * ".tex")
        end
            

        ens = load_ensemble(line)
        meas = nothing

        analysis = deepcopy(ens.global_metadata)
        contains_nans = ens.global_metadata[:nan_confs] != []
        
        if ensemble_data != nothing && "measurements" in keys(ensemble_data)
            meas = load_measurements(ensemble_data["measurements"], β = ens.global_metadata[:β], csw = ens.global_metadata[:csw])
            contains_nans = false
        end
	therm = get_key_or_nothing(ensemble_data, "therm")
        if !contains_nans
            if therm == nothing
                if ens.global_metadata[:nconfs] > 1100
                    therm = 1000
                elseif ens.global_metadata[:nconfs] > 500
                    therm = 400
                else
                    therm = 1
                end
            end
            measurements = ens
            if meas != nothing
                thermalise!(meas, therm)
                measurements = meas
            end

            thermalise!(ens, therm)
                
            analysis[:therm] = first(ens.analysis).confno
            corrs = [:g5_folded, :gk_folded, :id_folded]
            T = ens.global_metadata[:geometry][0]
            L = ens.global_metadata[:geometry][1]
            T_middle = T ÷ 2
            DEFAULT_FIT_WINDOW = (T_middle - NO_FIT_POINTS):T_middle
            function get_fit_window(dict, key)
                fit_window = get_key_or_nothing(dict, key)
                if fit_window == nothing
                    return DEFAULT_FIT_WINDOW
                end
                return fit_window
            end
            try
                @info "Plotting plaquette..."
                plot_plaquette(ens)
                save_figure("plaquette")
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "PCAC mass..."
                bs, tune = get_binsize_tune(ensemble_data, "pcac_binsize")
                fit_window = get_fit_window(ensemble_data, "pcac_fitwindow")
                analysis[:pcac] = bootstrap_effective_pcac(measurements.analysis, bs)
                plot_pcac_mass(measurements, bs)
                save_figure("pcac_mass")
                analysis[:m_pcac] = fit_pcac_mass(measurements, bs, fit_window)
                plot_pcac_fit(measurements, bs, fit_window, 1:T_middle)
                save_figure("pcac_mass_fit")
                if tune
                    tune_binsize_pcac_fit(measurements, DEFAULT_TUNE_BINSIZES, fit_window)
                    save_figure("pcac_mass_autocorrelations")
                end
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Fps..."
                bs, tune = get_binsize_tune(ensemble_data, "fps_binsize")
                fit_window = get_fit_window(ensemble_data, "fps_fitwindow")
                analysis[:effective_fps] = bootstrap_effective_fps(measurements.analysis, T, L, bs)
                plot_fps(measurements, bs)
                save_figure("fps")
                analysis[:fps] = fit_fps(measurements, bs, fit_window)
                plot_fps_fit(measurements, bs, fit_window, 1:T_middle)
                save_figure("fps_fit")
                if tune
                    tune_binsize_fps_fit(measurements, DEFAULT_TUNE_BINSIZES, fit_window)
                    save_figure("fps_autocorrelations")
                end
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Gps..."
                bs, tune = get_binsize_tune(ensemble_data, "gps_binsize")
                fit_window = get_fit_window(ensemble_data, "gps_fitwindow")
                analysis[:effective_gps] = bootstrap_effective_gps(measurements.analysis, ens.global_metadata[:geometry][1], bs)
                plot_gps(measurements, bs)
                save_figure("gps")
                analysis[:gps] = fit_gps(measurements, bs, fit_window)
                plot_gps_fit(measurements, bs, fit_window, 1:T_middle)
                save_figure("gps_fit")
                if tune
                    tune_binsize_gps_fit(measurements, DEFAULT_TUNE_BINSIZES, fit_window)
                    save_figure("gps_autocorrelations")
                end
            catch e
                @error "Failed!"
                @error e
            end
            for corr in corrs
                try
                    @info String(corr) * "..."
                    bs, tune = get_binsize_tune(ensemble_data, String(corr)*"_binsize")
                    fit_window = get_fit_window(ensemble_data, String(corr)*"_fitwindow")
                    analysis[Symbol("effective_"* String(corr))] = bootstrap_effective_mass(measurements.analysis, corr, bs)
                    plot_effective_mass(measurements, corr, bs)
                    save_figure("effective_mass_" * String(corr))
                    analysis[corr] = fit_effective_mass(measurements, corr, bs, fit_window)
                    plot_effective_mass_fit(measurements, corr, bs, fit_window, 1:T_middle)
                    save_figure("effective_mass_" * String(corr) * "_fit")
                    if tune
                        tune_effective_mass_fit(measurements, corr, DEFAULT_TUNE_BINSIZES, fit_window)
                        save_figure(String(corr)*"_autocorrelations")
                    end
                catch e
                    @error "Failed!"
                    @error e
                end
            end
            try
                @info "mv/mpi..."
                analysis[:effective_mass_ratio_mv_mpi] = bootstrap_effective_mass_ratio(measurements.analysis, :gk_folded, :g5_folded, bs)
                plot_effective_mass_ratio(measurements, :gk_folded, :g5_folded, bs)
                save_figure("effective_ratio_mv_mpi")
                analysis[:ratio_mv_mpi] = [last(analysis[:effective_mass_ratio_mv_mpi][1]), last(analysis[:effective_mass_ratio_mv_mpi][2])]
            catch e
                @error "Failed!"
                @error e
            end
            if wf != nothing
                @info "Wilson flow..."
                try
                    @info "t²E..."
                    plot_t2e(wf, binsize = bs)
                    save_figure("t2e")
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "W..."
                    plot_w(wf, binsize = bs)
                    save_figure("W")
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "Topological Charge..."
                    plot_tc(wf)
                    save_figure("tc")
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "Topological Charge Histogram..."
                    plot_tc_hist(wf)
                    save_figure("tc_hist")
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "Calculating w0..."
                    bs = get_binsize(ensemble_data, "w_binsize")
                    w0 = auto_w0(wf, binsize = bs)
                    analysis[:w0] = w0
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "Calculating mpiw0..."
                    analysis[:mpiw0] = propagate_product(analysis[:w0], analysis[:g5_folded])
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "Calculating mpiw0..."
                    analysis[:mpiw02] = propagate_product(analysis[:mpiw0], analysis[:mpiw0])
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "Calculating mpcacw0..."
                    analysis[:mpcacw0] = propagate_product(analysis[:w0], analysis[:m_pcac])
                catch e
                    @error "Failed!"
                    @error e
                end
            end
            try
                @info "Calculating mpi^2..."
                analysis[:mpi2] = propagate_product(analysis[:g5_folded], analysis[:g5_folded])
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Calculating mpiL..."
                analysis[:mpiL] = analysis[:g5_folded][1]*L
            catch e
                @error "Failed!"
                @error e
            end
        end
        YAML.write_file(ensemble_path * "analysis.yml", analysis)
        if contains_nans
            touch(ensemble_path * "CONTAINS_NANS")
        end
        return analysis
        
    end
    
end

for line in keys(list_of_ensembles)
    @info "Processing " * line
    analysis = process_ensemble(line, list_of_ensembles[line])
    ensembles[line] = analysis
end


YAML.write_file(OUTPUT_DIRECTORY * "ensembles.yml", ensembles)
