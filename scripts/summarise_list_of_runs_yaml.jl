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

#RIPREP_COMMIT = read(`git rev-parse --short HEAD`, String)

OUTPUT_DIRECTORY = "./" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS") * "/"
DEFAULT_BINSIZE = 10
NO_FIT_POINTS = 4
DEFAULT_TUNE_BINSIZES = collect(1:6:80)
TIKZ = false

NBOOT_LARGE = 10000
NBOOT_TUNING = 100

WF_REF = 1.0

ensembles = Dict()

ensure_directory_exists(OUTPUT_DIRECTORY)

list_of_ensembles = YAML.load_file(ARGS[1])

if TIKZ
    pgfplotsx()
else
    gr()
end

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
#    ensemble_path = OUTPUT_DIRECTORY * replace(line, "/" => "_")[2:end] * "/"
#    ensure_directory_exists(ensemble_path)
#    with_logger(FileLogger(ensemble_path * "analysis.log")) do
    with_logger(FileLogger(OUTPUT_DIRECTORY * "analysis.log")) do
        wf = nothing
        try
            wf = ensemble_data["wf"]
            wf = load_wilsonflow(wf)
            @info "Successfully loaded Wilson flow data"
        catch e
        end
        function save_figure(name)
            if TIKZ
                savefig(OUTPUT_DIRECTORY * name * ".tikz")

                s = read(open(OUTPUT_DIRECTORY * name * ".tikz"), String)
                open(name * ".tex", "w") do f
                    write(f, Plots.pgfx_preamble() * "\n \\begin{document} \n " * s * "\n \\end{document} \n")
                end
            else
                savefig(OUTPUT_DIRECTORY * name * ".pdf")
            end
            
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
            try
                @info "Plotting plaquette..."
                plot_plaquette(ens)
                save_figure("plaquette")
            catch e
                @error "Failed!"
                @error e
            end
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

            if therm < max(analysis[:nmpi_changes]...)
                thermalise!(ens, max(analysis[:nmpi_changes]...))
            else
                thermalise!(ens, therm)
            end

            try
                analysis[:mvm] = mean(ens.analysis.mvm)
            catch e
            end
            try
                analysis[:time] = mean(ens.analysis.time)
            catch e
            end
            
            thermalise!(ens, therm)

            try
                @info "Plotting plaquette thermalised..."
                plot_plaquette(ens)
                save_figure("plaquette_therm")
            catch e
                @error "Failed!"
                @error e
            end
                
            analysis[:therm] = first(ens.analysis).confno
            analysis[:acceptance] = mean(ens.analysis.accepted)
            corrs = [:g5_folded, :gk_folded, :id_folded]
            T = ens.global_metadata[:geometry][0]
            L = ens.global_metadata[:geometry][1]
            T_middle = T ÷ 2
            DEFAULT_FIT_WINDOW = (T_middle - NO_FIT_POINTS):(T_middle - 1)
            DEFAULT_PLOT_WINDOW = 3:(T_middle - 1)
            function get_fit_window(dict, key)
                fit_window = get_key_or_nothing(dict, key)
                if fit_window == nothing
                    return DEFAULT_FIT_WINDOW
                end
                return string_to_unitrange(fit_window)
            end
            function get_plot_window(dict, key)
                plot_window = get_key_or_nothing(dict, key)
                if plot_window == nothing
                    return DEFAULT_PLOT_WINDOW
                end
                return string_to_unitrange(plot_window)
            end
            try
                @info "Plotting fundamental Polyakov..."
                plot_fundamental_polyakov(ens)
                save_figure("fundamental_polyakov")
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Plotting fundamental Polyakov Histogram..."
                plot_fundamental_polyakov_hist(ens)
                save_figure("fundamental_polyakov_hist")
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Plotting adjoint Polyakov..."
                plot_adjoint_polyakov(ens)
                save_figure("adjoint_polyakov")
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Plotting adjoint Polyakov Histogram..."
                plot_adjoint_polyakov_hist(ens)
                save_figure("adjoint_polyakov_hist")
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Plotting maximum eigenvalue..."
                res = plot_maxeig(ens)
                if res != nothing
                    save_figure("maxeig")
                end
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "Plotting minimum eigenvalue..."
                res = plot_mineig(ens)
                if res != nothing
                    save_figure("mineig")
                end
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "PCAC mass..."
                bs, tune = get_binsize_tune(ensemble_data, "pcac_binsize")
                fit_window = get_fit_window(ensemble_data, "pcac_fitwindow")
                plot_window = get_plot_window(ensemble_data, "pcac_plotwindow")
                analysis[:pcac] = bootstrap_effective_pcac(measurements.analysis, bs)
                plot_pcac_mass(measurements, bs, plot_window)
                save_figure("pcac_mass")
                analysis[:m_pcac] = fit_pcac_mass(measurements, bs, fit_window)
                analysis[:pcac_binsize] = bs
                analysis[:pcac_fitwindow] = fit_window
                plot_pcac_fit(measurements, bs, fit_window, plot_window)
                save_figure("pcac_mass_fit")
                if tune
                    tune_binsize_pcac_fit(measurements, DEFAULT_TUNE_BINSIZES, fit_window, nboot = NBOOT_TUNING)
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
                plot_window = get_plot_window(ensemble_data, "fps_plotwindow")
                analysis[:effective_fps] = bootstrap_effective_fps(measurements.analysis, T, L, bs)
                plot_fps(measurements, bs, plot_window)
                save_figure("fps")
                analysis[:fps] = fit_fps(measurements, bs, fit_window)
                analysis[:fps_binsize] = bs
                analysis[:fps_fitwindow] = fit_window
                plot_fps_fit(measurements, bs, fit_window, plot_window)
                save_figure("fps_fit")
                if tune
                    tune_binsize_fps_fit(measurements, DEFAULT_TUNE_BINSIZES, fit_window, nboot = NBOOT_TUNING)
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
                plot_window = get_plot_window(ensemble_data, "gps_plotwindow")
                analysis[:effective_gps] = bootstrap_effective_gps(measurements.analysis, ens.global_metadata[:geometry][1], bs)
                plot_gps(measurements, bs, plot_window)
                save_figure("gps")
                analysis[:gps] = fit_gps(measurements, bs, fit_window)
                analysis[:gps_binsize] = bs
                analysis[:gps_fitwindow] = fit_window
                plot_gps_fit(measurements, bs, fit_window, plot_window)
                save_figure("gps_fit")
                if tune
                    tune_binsize_gps_fit(measurements, DEFAULT_TUNE_BINSIZES, fit_window, nboot = NBOOT_TUNING)
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
                    plot_window = get_plot_window(ensemble_data, String(corr)*"_plotwindow")
                    analysis[Symbol("effective_"* String(corr))] = bootstrap_effective_mass(measurements.analysis, corr, bs, range=plot_window)
                    plot_effective_mass(measurements, corr, bs, plot_window)
                    save_figure("effective_mass_" * String(corr))
                    analysis[corr] = fit_effective_mass(measurements, corr, bs, fit_window)
                    analysis[Symbol(String(corr) * "_binsize")] = bs
                    analysis[Symbol(String(corr) * "_fitwindow")] = bs
                    plot_effective_mass_fit(measurements, corr, bs, fit_window, plot_window)
                    save_figure("effective_mass_" * String(corr) * "_fit")
                    if tune
                        tune_effective_mass_fit(measurements, corr, DEFAULT_TUNE_BINSIZES, fit_window, nboot = NBOOT_TUNING)
                        save_figure(String(corr)*"_autocorrelations")
                    end
                catch e
                    @error "Failed!"
                    @error e
                end
            end
            try
                @info "mv/mpi naive..."
                bs, tune = get_binsize_tune(ensemble_data, "ratio_mv_mpi_binsize")
                fit_window = get_fit_window(ensemble_data, "ratio_mv_mpi_fitwindow")
                plot_window = get_plot_window(ensemble_data, "ratio_mv_mpi_plotwindow")
                analysis[:effective_mass_ratio_mv_mpi] = bootstrap_effective_mass_ratio(measurements.analysis, :gk_folded, :g5_folded, bs)
                analysis[:ratio_mv_mpi_naive] = fit_effective_mass_ratio_naive(measurements.analysis, :gk_folded, :g5_folded, bs, fit_window)
                plot_effective_mass_ratio_fit_naive(measurements, :gk_folded, :g5_folded, bs, fit_window, plot_window, nboot=NBOOT_LARGE)
                save_figure("effective_ratio_mv_mpi_fit")
                analysis[:ratio_mv_mpi_naive_binsize] = bs
                analysis[:ratio_mv_mpi_naive_fitwindow] = fit_window
                if tune
                    tune_effective_mass_ratio_fit_naive(measurements, :gk_folded, :g5_folded, DEFAULT_TUNE_BINSIZES, fit_window, nboot=NBOOT_TUNING)
                    save_figure("ratio_mv_mpi_autocorrelations")
                end
            catch e
                @error "Failed!"
                @error e
            end
            try
                @info "mv/mpi correlated..."
                bs_g5, tune_g5 = get_binsize_tune(ensemble_data, "g5_binsize")
                fit_window_g5 = get_fit_window(ensemble_data, "g5_fitwindow")
                bs_gk, tune_gk = get_binsize_tune(ensemble_data, "gk_binsize")
                fit_window_gk = get_fit_window(ensemble_data, "gk_fitwindow")
		analysis[:ratio_vs_mpi_correlated] = fit_ratio_correlated(measurements, :g5_folded, :gk_folded, bs_g5, fit_window_g5, bs_gk, fit_window_gk)
            catch e
                @error "Failed!"
            end
            
            if wf != nothing
                @info "Wilson flow..."
                try
                    @info "t²E..."
                    bs, tune = get_binsize_tune(ensemble_data, "t2e_binsize")
                    plot_t2e(wf, binsize = bs)
                    save_figure("t2e")
                catch e
                    @error "Failed!"
                    @error e
                end
                try
                    @info "W..."
                    bs, tune = get_binsize_tune(ensemble_data, "w_binsize")
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
                    bs, tune = get_binsize_tune(ensemble_data, "w0_binsize")
                    w0 = auto_w0(wf, binsize = bs)
                    analysis[:w0] = w0
                    analysis[:w0_binsize] = bs
		    analysis[:smearing_radius] = w0*sqrt(8)/L
                  if tune
                      tune_binsize_auto_w0(wf, DEFAULT_TUNE_BINSIZES, nboot = NBOOT_TUNING)
                      save_figure("w0_autocorrelations")
                  end
                catch e
                    @error "Failed!"
                    @error e
                end
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
    try
        analysis = process_ensemble(line, list_of_ensembles[line])
        ensembles[line] = analysis
    catch e
        @error "Error processing ensemble"
    end
end


YAML.write_file(OUTPUT_DIRECTORY * "ensembles.yml", ensembles)
