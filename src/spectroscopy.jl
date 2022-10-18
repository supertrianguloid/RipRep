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

function fit_cosh(ens, correlator, range; fold = :none, nstates = 1, p0 = :random)
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
        thermalise_bin!(ensemble, thermalisation, i);
        error = double_bootstrap_fits(ensemble, correlator, trange, nstates = nstates, p0 = p0, bs = bs)
        push!(errors, error)
    end
    errors = vcat(errors...)
    plot(binrange, errors[:, 2], yerr = errors[:, 4], title = _ensemble_to_latex_string(ensemble), label = "Error on the error in the mass for " * string(correlator))
end

function plot_mpcac(ens; nboot = 1000, folded = true, _bang = false)
    mpcac = pcac_mass(ens, folded = folded, nboot = nboot)
    μ = mean(mpcac)
    σ = std(mpcac)    

    plot_func = _bang ? plot! : plot

    plot_func(1:length(μ), μ, yerr=σ, label="\$m_{pcac}(\\tau)\$", title = _ensemble_to_latex_string(ens))
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
                push!(mpcac_inner, 0.5*mean(ens.analysis[bs2, :dg5_g0g5_re_folded])./mean(ens.analysis[bs2, :g5_folded])[1:end-1])
            else
                push!(mpcac_inner, 0.5*mean(ens.analysis[bs2, :dg5_g0g5_re])./mean(ens.analysis[bs2, :g5])[1:end-1])
            end
        end
        push!(mpcac_outer, mean(mpcac_inner))
    end
    return mpcac_outer
end

function pcac_mass_naive(ens::Ensemble; folded = true)
    
end

function pcac_mass(ens::Ensemble; folded = true, nboot = 1000, shift = 0)
    mpcac = []
    
    for n in 1:nboot
        bs = rand(1:nrow(ens.analysis), nrow(ens.analysis))
        if folded
            push!(mpcac, 0.5*mean(ens.analysis[bs, :dg5_g0g5_re_folded])./(mean(ens.analysis[bs, :g5_folded])[(1 + shift):end-1+shift]))
        else
            push!(mpcac, 0.5*mean(ens.analysis[bs, :dg5_g0g5_re])./(mean(ens.analysis[bs, :g5])[(1 + shift):end-1+shift]))
        end
    end
    return mpcac
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