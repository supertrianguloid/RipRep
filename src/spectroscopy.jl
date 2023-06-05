using Base: parentmodule_before_main
using Plots
using Statistics
using Roots
using LsqFit

include("parser.jl")
include("utilities.jl")

TITLE_FONT_SIZE = 13

function plot_plaquette(ens; range = :default, _bang = false)
    if range == :default
        range = ens.analysis[:, :confno]
    end
    plot_func = _bang ? plot! : plot
    plot_func(range, ens.analysis[in.(ens.analysis.confno, (range,)), :plaquette], label = "Plaquette", title = _ensemble_to_latex_string(ens), titlefontsize = TITLE_FONT_SIZE)
    xlabel!("Configuration")
end

function plot_plaquette!(ens; range = :default)
    plot_plaquette(ens, range=range, _bang=true)
end

function plot_correlator(ens, corr, range = :default; binsize = 1, log::Bool=true, _bang = false)
    if range == :default
        range = eachindex(ens.analysis[1, corr])
    end
    correlator = ens.analysis[:, corr]
    plot_func = _bang ? plot! : plot
    plot_func(Array(range), parent(mean(correlator)[range]), yerr = parent((std(correlator, corrected=true)/sqrt(length(correlator)))[range]), yaxis = log ? :log : :identity, label = String(corr), title = title = _ensemble_to_latex_string(ens), titlefontsize = TITLE_FONT_SIZE)
    xlabel!("\$τ\$")
end

function h(T, τ, E₁, E₂)
    return exp(-E₁*τ - E₂*(T - τ)) + exp(-E₂*τ - E₁*(T - τ))
end

function effective_pcac(df::DataFrame)
    T = length(df[1, :g5])
    meff = effective_mass(mean(df[:, :g5_folded]), T)
    return  0.5 .* meff./sinh.(meff) .* mean(df[:, :dg5_g0g5_re_folded])./(mean(df[:, :g5_folded])[1:end])
end

function bootstrap_effective_pcac(df::DataFrame, binsize; binmethod = :randomsample, nboot = 100)
    pcac = []
    for i in 1:nboot
        data = get_subsample(df, binsize, method=binmethod)
        push!(pcac, effective_pcac(data))
    end
    return mean(pcac), std(pcac)
end

function fit_pcac_mass(ens, binsize, fitting_range; binmethod = :randomsample, nboot = 100)
    fit = []
    for i in 1:nboot
        subsample = get_subsample(ens.analysis, binsize, method=binmethod)
        μ, σ = pcac_effective_mass(subsample, 1, binmethod=binmethod, nboot=nboot)
        push!(fit, only(fit_const(fitting_range, μ, σ).param))
    end
    μ = mean(fit)
    σ = std(fit)
    return μ, σ
end

function effective_mass(correlator, T)
    meffs = []
    eq(m, τ) = h(T, τ - 1, 0, m)/h(T, τ, 0, m) - correlator[τ - 1]/correlator[τ]
    for τ in eachindex(correlator)[1:end]
        push!(meffs, only(find_zeros(m -> eq(m, τ), 0, 100)))
    end
    return meffs               
end

function effective_gps(correlator, T, L)
    timeslices = eachindex(correlator)[1:end]
    meff = effective_mass(correlator, T)
    return [sqrt(meff[t]*correlator[t]/h(T, t, 0, meff[t])*L^3) for t in timeslices]
end

function bootstrap_effective_gps(df::DataFrame, L, binsize; binmethod = :randomsample, nboot = 1000)
    gps = []
    T = length(df[1, :g5])
    
    for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        c = mean(subsample[:, :g5_folded])
        push!(gps, effective_gps(c, T, L))
    end
    return mean(gps), std(gps)
end

function effective_fps(df::DataFrame, T, L)
    meff = effective_mass(mean(df[:, :g5_folded]), T)[1:end]
    gps = effective_gps(mean(df[:, :g5_folded]), T, L)[1:end]
    pcac = effective_pcac(df)
    return 2 .* (pcac./(meff.^2)).*gps
end

function bootstrap_effective_fps(df::DataFrame, T, L, binsize; binmethod = :randomsample, nboot = 1000)
    fps = []
    for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        push!(fps, effective_fps(subsample, T, L))
    end
    return mean(fps), std(fps)
end

function bootstrap_effective_mass(df::DataFrame, corr, binsize; binmethod = :randomsample, nboot = 1000)
    meff_boot = []
    
    T = length(df[1, :g5])
    
    for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        c = mean(subsample[:, corr])
        
        effective_masses = effective_mass(c, T)
        push!(meff_boot, effective_masses)
    end
    μ = mean(meff_boot)
    σ = std(meff_boot)
    return μ, σ
end

function fit_effective_mass(ens, corr, binsize, fitting_range; binmethod = :randomsample, range = :default, nboot = 1000)
    fit = []
    for i in 1:nboot
        subsample = get_subsample(ens.analysis, binsize, method=binmethod)
        μ, σ = bootstrap_effective_mass(subsample, corr, 1, binmethod=binmethod, range=range, nboot=nboot)
        push!(fit, only(fit_const(fitting_range, μ, σ).param))
               
    end
    μ = mean(fit)
    σ = std(fit)
    return μ, σ
end

#REWRITE
function plot_fit(ens::Ensemble, correlator::Symbol, trange; plotrange = :default, nstates = 1, log = true, p0 = :random)
    T = ens.global_metadata[:geometry][0]
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

#REWRITE
function bootstrap_fits(ens, correlator, trange; nstates = 1, bs = 100, p0 = :random)
    fits = []
    residuals = []
    stats = nrow(ens.analysis)
    T = ens.global_metadata[:geometry][0]
    for i ∈ 1:bs
        data = ens.analysis[rand(1:stats, stats), correlator]
        fit = _fit_helper(data, trange, nstates, T, p0)
        push!(fits, fit.param)
        push!(residuals, fit.resid)
    end
    return hcat(mean(fits)..., std(fits)...), hcat(residuals...)'
end
