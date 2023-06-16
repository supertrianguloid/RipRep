using YAML
include("spectroscopy.jl")

OUTPUT_DIRECTORY = "/home/laurence/projects/RipRep/output/"

function process_yaml_ensemble(path::String)
    yaml = YAML.load_file(path)
    overall_results = Dict()
    for ensemble in keys(yaml)
        ensemble_path = OUTPUT_DIRECTORY * replace(ensemble, "/" => "_")[2:end] * "/"
        ensure_directory_exists(ensemble_path)
        save_figure(name) = savefig(ensemble_path * name)
        save_results(results) = write(ensemble_path * "results.yml", YAML.write(results))
        results = Dict()
        metadata = yaml[ensemble]
        nboot = metadata["global_parameters"]["bootstrap_samples"]
        binmethod = metadata["global_parameters"]["bin_method"]
        ens = load_ensemble(ensemble)
        thermalise!(ens, metadata["global_parameters"]["thermalisation"])
        @info "Effective Masses..."
        for corr in keys(metadata["effective_masses"])
            params = metadata["effective_masses"][corr]
            @show corr
            plot_effective_mass_fit(ens, Symbol(corr), params["binsize"], string_to_unitrange(params["fitting_range"]), string_to_unitrange(params["plotting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
            save_figure("effective_mass_" * corr * ".pdf")
            results[corr] = fit_effective_mass(ens, corr, params["binsize"], string_to_unitrange(params["fitting_range"]), nboot=nboot, binmethod = Symbol(binmethod))
        end
        @info "PCAC..."
        if("pcac" in keys(metadata))
            params = metadata["pcac"]
            plot_pcac_fit(ens, params["binsize"], string_to_unitrange(params["fitting_range"]), string_to_unitrange(params["plotting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
            save_figure("pcac_mass.pdf")
            results[:pcac] = fit_pcac_mass(ens, params["binsize"], string_to_unitrange(params["fitting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
        end
        @info "Fps..."
        if("fps" in keys(metadata))
            params = metadata["fps"]
            plot_fps_fit(ens, params["binsize"], string_to_unitrange(params["fitting_range"]), string_to_unitrange(params["plotting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
            save_figure("fps.pdf")
            results[:fps] = fit_fps(ens, params["binsize"], string_to_unitrange(params["fitting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
        end
        @info "Gps..."
        if("gps" in keys(metadata))
            params = metadata["gps"]
            plot_gps_fit(ens, params["binsize"], string_to_unitrange(params["fitting_range"]), string_to_unitrange(params["plotting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
            save_figure("gps.pdf")
            results[:gps] = fit_fps(ens, params["binsize"], string_to_unitrange(params["fitting_range"]), binmethod = Symbol(binmethod), nboot = nboot)
        end
        save_results(results)
        overall_results[ensemble] = results
    end
    return overall_results
end

