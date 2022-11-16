using Plots
using Statistics
using LsqFit

include("parser.jl")
include("utilities.jl")

function plot_plaquette(ens; range = :default, _bang = false)
    if range == :default
        range = ens.analysis[:, :conf_no]
    end
    plot_func = _bang ? plot! : plot
    plot_func(range, ens.analysis[in.(ens.analysis.conf_no, (range,)), :plaquette], label="Plaquette", title = _ensemble_to_latex_string(ens))
    xlabel!("Configuration")
end

function plot_plaquette!(ens; range = :default)
    plot_plaquette(ens, range=range, _bang=true)
end

function default_thermalise_bin!(ens; bin = false)
    dtype = typeof(ens)
    if(dtype == Ensemble)
        N = only(ens.global_metadata.Confs)
    else
        N = nrow(ens.data)
    end
    if N > 1000
        therm = 800
        bs = 100
    elseif N > 750
        therm = 500
    else
        therm = 1
    end

    if N > 3000
        bs = 100
    elseif N > 1000
        bs = 10
    else
        bs = 5
    end
    if bin == false
        bs = 1
    end
    thermalise_bin!(ens, therm, bs)
end

function _cosh_model(nstates, T, τ, params)
    if nstates == 1
        A, m = params
        return @. A * cosh(m * (τ - (T/2)))
    elseif nstates == 2
        A₁, m₁, A₂, m₂ = params
        return @. A₁ * cosh(m₁ * (τ - (T/2))) + A₂ * cosh(m₂ * (τ - (T/2)))
    end
end

function plot_correlator(ens, corr, range = :default; log::Bool=true, _bang = false)
    if range == :default
        range = eachindex(ens.analysis[1, corr])
    end
    correlator = ens.analysis[:, corr]
    plot_func = _bang ? plot! : plot
    plot_func(Array(range), parent(mean(correlator)[range]), yerr = parent((std(correlator, corrected=true)/sqrt(length(correlator)))[range]), yaxis = log ? :log : :identity, label = String(corr), title = title = _ensemble_to_latex_string(ens))
    xlabel!("\$τ\$")
end

function plot_correlator!(ens, corr, range = :default; log::Bool=true)
    plot_correlator(ens, corr, range, log = log, _bang = true)
end

function fit_cosh(ens, correlator, range; nstates = 1, p0 = :random)
    T = ens.global_metadata.T
    data = ens.analysis[:, correlator]
    return _fit_helper(data, range, nstates, T, p0)
end

function _fit_helper(data, trange, nstates, T, p0)
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    μ = mean(data)
    w = 1 ./ var(data)
    if(p0 == :random)
        p0 = abs.(randn(nstates * 2)) .+ 1
    end
    lb = zeros(nstates * 2)
    fit = curve_fit(model, trange, μ[trange], w[trange], p0, lower=lb)
    if(findmax(fit.resid.^2)[1] > 1)
        #@error "Bad fit."
    end
    return fit
end

function plot_residuals(ens::Ensemble, correlator::Symbol, tmax ;nstates = 1, p0 = :random)
    resid = []
    for i in 1:tmax-3
        fit = fit_cosh(ens, correlator, i:tmax, nstates = nstates, p0 = p0)
        push!(resid, sum((fit.resid).^2))
    end
    plot(1:tmax-3, resid, yaxis = :log)
end



function plot_fit(ens::Ensemble, correlator::Symbol, trange; plotrange = :default, nstates = 1, log = true, p0 = :random)
    T = ens.global_metadata.T
    fit = fit_cosh(ens, correlator, trange, nstates = nstates, p0 = p0)
    println("Residuals: ", sum(fit.resid.^2))
    println("Fit parameters: ", fit.param)
    corr = ens.analysis[:, correlator]
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    if(plotrange == :default)
        plotrange = 0:length(mean(corr)) - 1
    end
    #plot(plotrange, parent(mean(corr))[plotrange .+ 1], yerr = parent(std(corr, corrected=true)/sqrt(length(corr)))[plotrange .+ 1], yaxis = log ? :log : :identity, label = String(correlator), title = _ensemble_to_latex_string(ens))
    #xlabel!("Imaginary Time \$τ\$")
    plot_correlator(ens, correlator, plotrange, log=log)
    #annotate!((0.25,0.25), string(fit.param))
    plot!(trange, model(trange, fit.param), label = "$nstates state fit")    
end

function bootstrap_fits(ens, correlator, trange; nstates = 1, bs = 100, p0 = :random)
    fits = []
    residuals = []
    stats = nrow(ens.analysis)
    T = ens.global_metadata.T
    for i ∈ 1:bs
        data = ens.analysis[rand(1:stats, stats), correlator]
        fit = _fit_helper(data, trange, nstates, T, p0)
        push!(fits, fit.param)
        push!(residuals, fit.resid)
    end
    return hcat(mean(fits)..., std(fits)...), hcat(residuals...)'
end

function double_bootstrap_fits(ens, correlator, trange; nstates = 1, bs = 100, p0 = :random)
    errors = []
    stats = nrow(ens.analysis)
    T = ens.global_metadata.T
    for i ∈ 1:bs
        data = ens.analysis[rand(1:stats, stats), correlator]
        fits_inner = []
        residuals_inner = []
        for j in 1:bs
            data2 = data[rand(1:stats, stats)]
            fit = _fit_helper(data2, trange, nstates, T, p0)
            push!(fits_inner, fit.param)
            push!(residuals_inner, fit.resid)
        end
        push!(errors, std(fits_inner))
    end
    return hcat(mean(errors)..., std(errors)...)
end

function plot_error(ensemble, correlator, trange, thermalisation, binrange; nstates = 1, bs = 20, p0 = :random, method = :equal, _bang = false)
    plot_func = _bang ? plot! : plot
    errors = []
    for i in binrange
        thermalise_bin!(ensemble, thermalisation, i, method = method);
        error = double_bootstrap_fits(ensemble, correlator, trange, nstates = nstates, p0 = p0, bs = bs)
        push!(errors, error)
    end
    errors = vcat(errors...)
    plot_func(binrange, errors[:, 2], yerr = errors[:, 4], title = _ensemble_to_latex_string(ensemble), label = "Error on the error in the mass for " * string(correlator))
end

function plot_error!(ensemble, correlator, trange, thermalisation, binrange; nstates = 1, bs = 20, p0 = :random, method = :equal, _bang = true)
    plot_error(ensemble, correlator, trange, thermalisation, binrange; nstates = nstates, bs = bs, p0 = p0, method = method, _bang = true)
end

function plot_mpcac(ens; nboot = 1000, folded = true, _bang = false)
    mpcac = pcac_mass(ens.analysis, folded = folded, nboot = nboot)
    μ = mean(mpcac)
    σ = std(mpcac)    

    plot_func = _bang ? plot! : plot

    plot_func(1:length(μ), μ, yerr=σ, label="\$m_{pcac}(\\tau)\$", title = _ensemble_to_latex_string(ens))
end

function plot_pcac_fit(ens, fitrange; nboot = 1000, folded = true, _bang = false)
    fit = _pcac_fit(ens.analysis, fitrange, nboot = nboot, folded = folded)
    m = only(fit.param)
    println("Residuals: ", sum(fit.resid.^2))
    println("Fit parameters: ", fit.param)
    plot(fitrange, repeat([m],length(fitrange)))
    plot_mpcac!(ens, nboot = nboot, folded = folded)
end

function _pcac_fit(data, fitrange; nboot, folded = true)
    mpcac = pcac_mass(data, folded = folded, nboot = nboot)
    μ = mean(mpcac)
    w = 1 ./ var(mpcac)
    const_model(t, p) = p[1] .+ 0 .* t
    fit = curve_fit(const_model, Vector(fitrange), μ[fitrange], w[fitrange], [1.0])
    return fit
end

function pcac_fit_bootstrap(ens, fitrange; folded = true, nboot = 100)
    mpcac_outer = []
    for n in 1:nboot
        bs = rand(1:nrow(ens.analysis), nrow(ens.analysis))
        push!(mpcac_outer, only(_pcac_fit(ens.analysis[bs, :], fitrange, nboot = nboot, folded = folded).param))
    end
    return mean(mpcac_outer), std(mpcac_outer)
end

function plot_mpcac!(ens; nboot = 1000, folded = true)
    plot_mpcac(ens, nboot = nboot, folded = folded, _bang = true)
end

function pcac_mass_double_bootstrap(ens::Ensemble; folded = true, nboot = 1000)
    mpcac_outer = []

    for n in 1:nboot
        bs1 = rand(1:nrow(ens.analysis), nrow(ens.analysis))
        mpcac_inner = []
        for m in 1:nboot
            bs2 = rand(bs1, nrow(ens.analysis))
            if folded
                mass = 0.5*mean(ens.analysis[bs2, :dg5_g0g5_re_folded])./mean(ens.analysis[bs2, :g5_folded])[1:end-1]
                push!(mpcac_inner, mass)
            else
                push!(mpcac_inner, 0.5*mean(ens.analysis[bs2, :dg5_g0g5_re])./mean(ens.analysis[bs2, :g5])[1:end-1])
            end
        end
        push!(mpcac_outer, mean(mpcac_inner))
    end
    return mpcac_outer
end
function pcac_mass(data; folded = true, nboot = 1000, shift = 0)
    mpcac = []
    
    for n in 1:nboot
        bs = rand(1:nrow(data), nrow(data))
        push!(mpcac, pcac_values(data[bs, :], folded = folded, shift = shift))
    end
    return mpcac
end

function pcac_values(data; folded = true, shift = 0)
    if folded
        return 0.5*mean(data[:, :dg5_g0g5_re_folded])./(mean(data[:, :g5_folded])[(1 + shift):end-1+shift])
    else
        return 0.5*mean(data[:, :dg5_g0g5_re])./(mean(data[:, :g5])[(1 + shift):end-1+shift])
    end
end

# TODO: BROKEN
function plot_pcac_error(ensemble, fitrange, thermalisation, binrange; nboot = 20, folded = true, method = :equal)
    errors = []
    for i in binrange
        thermalise_bin!(ensemble, thermalisation, i, method = method);
        error = pcac_fit_bootstrap(ensemble, fitrange; folded, nboot)[2]
        push!(errors, error)
    end
    μ = mean(errors)
    σ = std(errors)
    plot(binrange, μ, yerr = σ, title = _ensemble_to_latex_string(ensemble), label = "Error on the error in the PCAC mass")
end

function plot_effective_mass(ens, corr, range = :default; nboot = 1000, log = true)
    meff_boot = []
    if range == :default
        range = 0:length(mean(ens.analysis[1,corr]))-1
    end
    for i in 1:nboot
        bs = rand(1:nrow(ens.analysis), nrow(ens.analysis))
        c = mean(ens.analysis[bs, corr])
        
        meffs = []
        T = only(ens.global_metadata.T)
        for τ in Array(eachindex(c))[3:end]
            push!(meffs, only(find_zeros(m -> (c[τ-1]/c[τ])*cosh(m * (τ - (T/2))) - cosh(m * ((τ - 1) - (T/2))), 0.0,5.0)))
        end
        push!(meff_boot, meffs)
    end
    μ = mean(meff_boot)
    σ = std(meff_boot)
    plot(eachindex(μ), μ, yerr = σ, yaxis = log ? :log : :identity, label = String(corr) * " effective mass", title = title = _ensemble_to_latex_string(ens))
end

function _fps(data, T; folded = true, piwindow, piinitial, pcacwindow, nboot = 100)
    vals = []
    correlator = folded ? :g5_folded : :g5
    for i in 1:nboot
        bs = rand(1:nrow(data), nrow(data))
        fit = _fit_helper(data[bs, correlator], piwindow, 1, T, piinitial)
        A, m = fit.param
        mpcac = only(_pcac_fit(data[bs, :], pcacwindow, nboot = nboot, folded = folded).param)
        push!(vals, 2*mpcac*sqrt(A*m)/(m^2))
    end
    return mean(vals), std(vals)
end

function fps(ens; folded = true, piwindow, piinitial, pcacwindow, nboot = 100)
    return _fps(ens.analysis, ens.global_metadata.T, folded = folded, nboot = nboot, piwindow = piwindow, piinitial = piinitial, pcacwindow = pcacwindow)
end

#TODO: Broken
function plot_fps_error(ensemble, therm, maxbins; folded = true, piwindow, piinitial, pcacwindow, nboot = 100, method = :equal)
    vals = []
    for i in 1:maxbins
        thermalise_bin!(ensemble, therm, i, method = method);
        push!(vals, fps(ensemble, folded = folded, nboot = nboot, piwindow = piwindow, piinitial = piinitial, pcacwindow = pcacwindow)[2])
    end
end