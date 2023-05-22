using Base: parentmodule_before_main
using Plots
using Statistics
using Roots

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

function default_thermalise_bin!(ens; bin = false)
    dtype = typeof(ens)
    if(dtype == Ensemble)
        N = ens.global_metadata[:nconfs]
    else
        N = nrow(ens.data)
    end
    if N > 1000
        therm = 800
        bs = 100
    elseif N > 750
        therm = 500
    else
        @warn "Short run (N < 750), not thermalising!"
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

function plot_correlator(ens, corr, range = :default; log::Bool=true, _bang = false)
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

function effective_mass(correlator, T)
    meffs = []
    eq(m, τ) = (correlator[τ]/correlator[τ + 1])*(cosh(m*(τ + 1 - T÷2))) - (cosh(m*(τ - T÷2)))
    for τ in eachindex(correlator)[1:end-1]
        push!(meffs, only(find_zeros(m -> eq(m, τ), 0, 100)))
    end
    return meffs               
end

function effective_mass_2(correlator, T)
    meffs = []
    eq(m, τ) = h(T, τ - 1, 0, m)/h(T, τ, 0, m) - correlator[τ - 1]/correlator[τ]
    for τ in eachindex(correlator)[1:end]
        push!(meffs, only(find_zeros(m -> eq(m, τ), 0, 100)))
    end
    return meffs               
end
    
function plot_effective_mass(ens, corr; range = :default, nboot = 1000, log = true, _bang = false)
    plot_func = _bang ? plot! : plot
    meff_boot = []

    if range == :default
        range = eachindex(ens.analysis[1, corr])[1:end-1]
    end
    for i in 1:nboot
        bs = rand(1:nrow(ens.analysis), nrow(ens.analysis))
        c = mean(ens.analysis[bs, corr])
        
        T = only(ens.global_metadata[:geometry][0])
        effective_masses = effective_mass(c, T)
        push!(meff_boot, effective_masses)
    end
    μ = mean(meff_boot)
    σ = std(meff_boot)
    @show μ
    @show σ
    plot_func(range, μ[range], yerr = σ[range], yaxis = log ? :log : :identity, title = String(corr) * " effective mass", legend = false)
end

function plot_effective_mass!(ens, corr; range = :default, nboot = 1000, log = true)
    plot_effective_mass(ens, corr, range=range, nboot=nboot, log=log, _bang=true)
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
