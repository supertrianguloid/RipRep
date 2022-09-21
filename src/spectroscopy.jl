using Plots
using Statistics
using LsqFit

theme(:dracula)
default(size = (900, 700))

include("parser.jl")
include("utilities.jl")

@info "RipRep spectroscopy code activated."


function plot_plaquette(ens, range = :default)
    if range == :default
        range = ens.analysis[:, :conf_no]
    end
    plot(range, ens.analysis[in.(ens.analysis.conf_no, (range,)), :plaquette], label="Plaquette", title = _ensemble_to_latex_string(ens))
    xlabel!("Configuration")
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

function plot_correlator(ens, corr; log::Bool=true)
    correlator = ens.analysis[:, corr]
    plot_correlator(ens, corr, 0:length(mean(correlator)) - 1, log=log)
end

function plot_correlator(ens, corr, range; log::Bool=true)
    correlator = ens.analysis[:, corr]
    plot(range, parent(mean(correlator)[range]), yerr = parent((std(correlator, corrected=true)/sqrt(length(correlator)))[range]), yaxis = log ? :log : :identity, label = String(corr), title = title = _ensemble_to_latex_string(ens))
    xlabel!("\$τ\$")
end

function fit_cosh(ens, correlator, range; fold = :none, nstates = 1, p0 = :random)
    T = ens.global_metadata.T
    data = ens.analysis[:, correlator]
    return _fit_helper(data, range, nstates, T, p0)
end

function _fit_helper(data, trange, nstates, T, p0)
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    μ = mean(data)
    w = 1 ./ var(data, corrected=true)
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

function plot_fit(ens::Ensemble, correlator::Symbol, trange; plotrange = :default, fold=:none, nstates = 1, log = true, p0 = :random)
    T = ens.global_metadata.T
    fit = fit_cosh(ens, correlator, trange, fold = fold, nstates = nstates, p0 = p0)
    println("Residuals: ", fit.resid)
    println("Fit parameters: ", fit.param)
    corr = ens.analysis[:, correlator]
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    if(plotrange == :default)
        plotrange = 0:length(mean(corr)) - 1
    end
    plot(plotrange, parent(mean(corr))[plotrange .+ 1], yerr = parent(std(corr, corrected=true)/sqrt(length(corr)))[plotrange .+ 1], yaxis = log ? :log : :identity, label = String(correlator), title = _ensemble_to_latex_string(ens))
    xlabel!("Imaginary Time \$τ\$")
    #annotate!((0.25,0.25), string(fit.param))
    plot!(trange, model(trange, fit.param), label = "$nstates state fit")    
end

function thermalise!(ens::Ensemble, ntherm)
    # TODO: Here would be a good place to add a check to make sure we got rid of all the missing data, and throw an error if not.
    ens.analysis = ens.data[ntherm:end, :]
end

function bin!(ens, binsize, method = :equal)
    data = ens.analysis
    if binsize == 1
        return
    end
    offset = nrow(data) % binsize
    if(offset != 0)
        data = data[offset:end, :]
    end
    newlen = nrow(data) ÷ binsize
    if method == :equal
        ens.analysis = data[[1 + (i - 1)*binsize for i in 1:newlen], :]
    end
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
            data = data[rand(1:stats, stats)]
            fit = _fit_helper(data, trange, nstates, T, p0)
            push!(fits_inner, fit.param)
            push!(residuals_inner, fit.resid)
        end
        push!(errors, std(fits_inner))
    end
    return hcat(mean(errors)..., std(errors)...)
end

function plot_error(ensemble, correlator, trange, thermalisation, binrange; nstates = 1, bs = 20, p0 = :random)
    errors = []
    for i in binrange
        thermalise!(ensemble, thermalisation);
        bin!(ensemble, i);
        error = double_bootstrap_fits(ensemble, correlator, trange, nstates = nstates, p0 = p0, bs = bs)
        push!(errors, error)
    end
    errors = vcat(errors...)
    plot(binrange, errors[:, 2], yerr = errors[:, 4], title = _ensemble_to_latex_string(ensemble), label = "Error on the error in the mass for " * string(correlator))
end