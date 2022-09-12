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
        range = ens.data[:, :conf_no]
    end
    plot(range, ens.data[in.(ens.data.conf_no, (range,)), :plaquette], label="Plaquette", title = _ensemble_to_latex_string(ens))
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

function plot_correlator(ens, corr; _log::Bool=true)
    correlator = ens.data[:, corr]
    plot_correlator(ens, corr, 0:length(mean(correlator)) - 1, _log=_log)
end

function plot_correlator(ens, corr, trange; _log::Bool=true)
    correlator = ens.data[:, corr]
    plot(trange, parent(mean(correlator)[trange]), yerr = parent((std(correlator, corrected=true)/sqrt(length(correlator)))[trange]), yaxis = _log ? :log : :identity, label = String(corr), title = title = _ensemble_to_latex_string(ens))
    xlabel!("\$τ\$")
end

function fit_cosh(ens, correlator, trange, nstates = 1, p0 = :random)
    T = ens.global_metadata.T
    data = ens.data[:, correlator]
    return _fit_helper(data, trange, nstates, T, p0)
end

function _fit_helper(data, trange, nstates, T, p0)
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    μ = mean(data)
    w = 1 ./ var(data, corrected=true)
    if(p0 == :random)
        p0 = abs.(randn(nstates * 2)) .+ 1
    end
    lb = zeros(nstates * 2)
    return curve_fit(model, trange, μ[trange], w[trange], p0, lower=lb)
end

function plot_fit(ens, correlator, trange, nstates = 1, log = true, p0 = :random)
    T = ens.global_metadata.T
    fit = fit_cosh(ens, correlator, trange, 1, p0)
    println(fit.resid)
    corr = ens.data[:, correlator]
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    plot(0:length(mean(corr)) - 1, parent(mean(corr)), yerr = parent(std(corr, corrected=true)/sqrt(length(corr))), yaxis = log ? :log : :identity, label = String(correlator), title = _ensemble_to_latex_string(ens))
    xlabel!("Imaginary Time \$τ\$")
    plot!(trange, model(trange, fit.param), label = "$nstates state fit")    
end