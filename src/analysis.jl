using DrWatson
using Plots
using Statistics
using LsqFit

theme(:dracula)

@quickactivate "RipRep"
include(srcdir("parser.jl"))

function acceptance(data)
    return sum(data[:, :accepted])/nrow(data)
end

function thermalise(data, ntherm)
    return data[ntherm:end, :]
end

function bin(data, binsize, method)
    offset = nrow(data) % binsize
    data = thermalise(data, offset)
    newlen = nrow(data) ÷ binsize
    if method == :equal
        return data[[1 + (i - 1)*binsize for i in 1:newlen], :]
    end
end

function plot_plaquette(ens)
    plot(1:nrow(ens.data), ens.data[:, :plaquette], label="Plaquette", title = only(ens.global_metadata.path))
    xlabel!("Conf #")
end
function plot_plaquette(ens, range)
    plot(range, ens.data[range, :plaquette], label="Plaquette", title = only(ens.global_metadata.path))
    xlabel!("Conf #")
end

function plot_correlator(ens, corr, log=true)
    correlator = ens.data[:, corr]
    plot(0:length(mean(correlator)) - 1, parent(mean(correlator)), yerr = parent(std(correlator, corrected=true)/sqrt(length(correlator))), yaxis = log ? :log : :identity, label = String(corr), title = only(ens.global_metadata.path))
    xlabel!("\$τ\$")
end

function plot_correlator(ens, corr, trange, log=true)
    correlator = ens.data[:, corr]
    plot(trange, parent(mean(correlator))[trange], yerr = parent(std(correlator, corrected=true)/sqrt(length(correlator)))[trange], yaxis = log ? :log : :identity, label = String(corr), title = only(ens.global_metadata.path))
    xlabel!("\$τ\$")
end

function _fit_helper(data, trange, nstates, T)
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    μ = mean(data)
    w = 1 ./ var(data, corrected=true)
    p0 = abs.(randn(nstates * 2)) .+ 1
    return curve_fit(model, trange, μ[trange], w[trange], p0)
end

function fit_cosh(ens, correlator, trange, nstates = 1)
    T = ens.global_metadata.T
    data = ens.data[:, correlator]
    return _fit_helper(data, trange, nstates, T)
end

function plot_fit(ens, correlator, trange, nstates = 1, log = true)
    T = ens.global_metadata.T
    fit = fit_cosh(ens, correlator, trange, 1)
    corr = ens.data[:, correlator]
    model(τ, params) = _cosh_model(nstates, T, τ, params)
    plot(0:length(mean(corr)) - 1, parent(mean(corr)), yerr = parent(std(corr, corrected=true)/sqrt(length(corr))), yaxis = log ? :log : :identity, label = String(correlator))
    xlabel!("\$τ\$")
    plot!(trange, model(trange, fit.param), label = "$nstates state fit")
    
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

function bootstrap_fits(ens, correlator, trange, nstates = 1, bs = 100)
    fits = []
    T = ens.global_metadata.T
    for i ∈ 1:bs
        data = ens.data[rand(1:bs, bs), correlator]
        push!(fits, _fit_helper(data, trange, nstates, T).param)
    end
    return hcat(mean(fits)..., std(fits)...)
end

function fits(ens, correlator, tmin, tmax, nstates = 1, bs = 100)
    fits = []
    for t in tmin:(tmax - 2*nstates)
        push!(fits, bootstrap_fits(ens, correlator, t:tmax, nstates, bs))
    end
    fits = vcat(fits...)
    plot(tmin:tmax - nstates*2, fits[:, 2], yerr = fits[:, 4], label = string(correlator) * " mass, \$\\tau_{max}\$ = $tmax", title = only(ens.global_metadata.path))
    xlabel!("Lower fitting range")
end