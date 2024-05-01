using Plots
using Statistics
using Roots
using LsqFit
using DataFrames
using Base.Threads: @spawn
using LaTeXStrings

include("parser.jl")
include("utilities.jl")

TITLE_FONT_SIZE = 13
NBOOT_DEFAULT = 20

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

function plot_correlator(ens, corr, range = :default; binsize = 1, nboot=NBOOT_DEFAULT, log::Bool=true, _bang = false)
    if range == :default
        range = eachindex(ens.analysis[1, corr])
    end
    correlator = ens.analysis[:, corr]
    plot_func = _bang ? plot! : plot
    plot_func(Array(range), parent(mean(correlator)[range]), yerr = parent(standard_error(correlator, binsize=binsize, nboot=nboot)[range]), yaxis = log ? :log : :identity, label = String(corr), title = title = _ensemble_to_latex_string(ens), titlefontsize = TITLE_FONT_SIZE)
    xlabel!("\$τ\$")
end

function h(T, τ, E₁, E₂)
    return exp(-E₁*τ - E₂*(T - τ)) + exp(-E₂*τ - E₁*(T - τ))
end

function effective_pcac(df::DataFrame)
    t = length(df[1, :g5])
    meff = effective_mass(mean(df[:, :g5_folded]), t)
    return  0.5 .* meff[1:end]./sinh.(meff[1:end]) .* mean(df[:, :dg5_g0g5_re_folded])./(mean(df[:, :g5_folded])[1:end])
end

function bootstrap_effective_pcac(df::DataFrame, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn try
                    effective_pcac(get_subsample(df, binsize, method=binmethod))
                catch error
                  missing
                end for i in 1:nboot]
    

    pcac = fetch.(res)
    nfailures = sum(ismissing.(pcac))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    pcac = [i for i in pcac if !ismissing(i)]
    return [mean(pcac), std(pcac)]
end

function bootstrap_effective_pcac_distribution(df::DataFrame, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn effective_pcac(get_subsample(df, binsize, method=binmethod)) for i in 1:nboot]
    pcac = fetch.(res)
    return pcac
end

function fit_pcac_mass(analysis::DataFrame, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn   try
                        subsample = get_subsample(analysis, binsize, method=binmethod)
                        μ, σ = bootstrap_effective_pcac(subsample, 1, binmethod=binmethod, nboot=nboot)
                        only(fit_const(fitting_range, μ, σ).param)
                    catch error
                        missing
                    end for i in 1:nboot]
    fit = fetch.(res)
    nfailures = sum(ismissing.(fit))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    fit = [i for i in fit if !ismissing(i)]

    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end
function fit_pcac_mass(ens::Ensemble, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    return fit_pcac_mass(ens.analysis, binsize, fitting_range, binmethod = binmethod, nboot = nboot)
end

function tune_binsize_pcac_fit(ens, binsizes, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(ens.analysis, bs)
            push!(fit_error, fit_pcac_mass(subsample, 1, fitting_range, binmethod = binmethod, nboot = nboot)[2])
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end

function effective_mass(correlator, T)
    meffs = Float64[]
    eq(m, τ) = h(T, τ - 1, 0, m)/h(T, τ, 0, m) - correlator[τ - 1]/correlator[τ]
    for τ in (first(eachindex(correlator))+1):last(eachindex(correlator))
        push!(meffs, only(find_zeros(m -> eq(m, τ), 0, 100)))
    end
    return meffs               
end

function effective_mass_ratio(corr_numerator, corr_denominator, T)
    return effective_mass(corr_numerator, T)./effective_mass(corr_denominator, T)
end


function effective_gps(correlator, T, L)
    timeslices = eachindex(correlator)[1:end]
    meff = effective_mass(correlator, T)
    return [sqrt(meff[t]*correlator[t]/h(T, t, 0, meff[t])*L^3) for t in timeslices]
end

function bootstrap_effective_gps(df::DataFrame, L, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    T = length(df[1, :g5])
    res = [@spawn try
                subsample = get_subsample(df, binsize, method=binmethod)
                c = mean(subsample[:, :g5_folded])
                effective_gps(c, T, L)
           catch e
               missing
           end for i in 1:nboot]
    gps = fetch.(res)
    nfailures = sum(ismissing.(gps))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    gps = [i for i in gps if !ismissing(i)]
    return [mean(gps), std(gps)]
end

function fit_gps(analysis::DataFrame, L, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn try
                subsample = get_subsample(analysis, binsize, method=binmethod)
                μ, σ = bootstrap_effective_gps(subsample, L, 1, binmethod=binmethod, nboot=nboot)
                only(fit_const(fitting_range, μ, σ).param)
           catch e
               missing
           end for i in 1:nboot]
    fit = fetch.(res)
    nfailures = sum(ismissing.(fit))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    fit = [i for i in fit if !ismissing(i)]
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function fit_gps(ens::Ensemble, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    return fit_gps(ens.analysis, L, binsize, fitting_range, binmethod = binmethod, nboot = nboot)
end

function tune_binsize_gps_fit(ens, binsizes, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(ens.analysis, bs)
            push!(fit_error, fit_gps(subsample, L, 1, fitting_range, binmethod = binmethod, nboot = nboot)[2])
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end

function effective_fps(df::DataFrame, T, L)
    meff = effective_mass(mean(df[:, :g5_folded]), T)[1:end]
    gps = effective_gps(mean(df[:, :g5_folded]), T, L)[1:end]
    pcac = effective_pcac(df)
    return 2 .* (pcac./(meff.^2)).*gps
end

function bootstrap_effective_fps(df::DataFrame, T, L, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn try
                    subsample = get_subsample(df, binsize, method=binmethod)
                    effective_fps(subsample, T, L)
           catch e
               missing
               end for i in 1:nboot]
    fps = fetch.(res)
    nfailures = sum(ismissing.(fps))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    fps = [i for i in fps if !ismissing(i)]
    return [mean(fps), std(fps)]
end

function fit_fps(ens::Ensemble, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    T = first(ens.global_metadata[:geometry])
    return fit_fps(ens.analysis, L, T, binsize, fitting_range, binmethod = binmethod, nboot = nboot)
end
function fit_fps(analysis::DataFrame, L, T, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn try
                    subsample = get_subsample(analysis, binsize, method=binmethod)
                    μ, σ = bootstrap_effective_fps(subsample, T, L, 1, binmethod=binmethod, nboot=nboot)
                    only(fit_const(fitting_range, μ, σ).param)
           catch e
               missing
               end for i in 1:nboot]
    fit = fetch.(res)
    nfailures = sum(ismissing.(fit))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    fit = [i for i in fit if !ismissing(i)]
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function tune_binsize_fps_fit(ens, binsizes, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    T = first(ens.global_metadata[:geometry])
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(ens.analysis, bs)
            push!(fit_error, fit_fps(subsample, L, T, 1, fitting_range, binmethod = binmethod, nboot = nboot)[2])
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end

function bootstrap_effective_mass(df::DataFrame, corr, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT, range=:all)
    
    T = length(df[1, :g5])
    
    res = [@spawn try
                    subsample = get_subsample(df, binsize, method=binmethod)
                    c = mean(subsample[:, corr])
                    if range != :all
                        c = c[(first(range)-1):last(range)]
                    end
                    effective_mass(c, T)
           catch e
               missing
               end for i in 1:nboot]
    meff_boot = fetch.(res)
    nfailures = sum(ismissing.(meff_boot))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    meff_boot = [i for i in meff_boot if !ismissing(i)]
    μ = mean(meff_boot)
    σ = std(meff_boot)
    return [μ, σ]
end

function bootstrap_effective_mass_ratio(df::DataFrame, corr_numerator, corr_denominator, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT, range=:all)
    
    T = length(df[1, :g5])
    
    res = [@spawn try
                    subsample = get_subsample(df, binsize, method=binmethod)
                    numerator = mean(subsample[:, corr_numerator])
                    denominator = mean(subsample[:, corr_denominator])
                    if range != :all
                        numerator = numerator[range]
                        denominator = denominator[range]
                    end

                    effective_mass_ratio(numerator, denominator, T)
           catch e
               missing
               end for i in 1:nboot]
    meff_boot = fetch.(res)
    nfailures = sum(ismissing.(meff_boot))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    meff_boot = [i for i in meff_boot if !ismissing(i)]
    μ = mean(meff_boot)
    σ = std(meff_boot)
    return [μ, σ]
end

function fit_effective_mass(analysis::DataFrame, corr, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn try
                      subsample = get_subsample(analysis, binsize, method=binmethod)
                      μ, σ = bootstrap_effective_mass(subsample, corr, 1, binmethod=binmethod, nboot=nboot, range=fitting_range)
                      only(fit_const(eachindex(μ), μ, σ).param)
           catch e
               missing
               end for i in 1:nboot]
    fit = fetch.(res)
    nfailures = sum(ismissing.(fit))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    fit = [i for i in fit if !ismissing(i)]
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function fit_effective_mass(ens::Ensemble, corr, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    return fit_effective_mass(ens.analysis, corr, binsize, fitting_range, binmethod = binmethod, nboot = nboot)
end

function tune_effective_mass_fit(ens, corr, binsizes, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(ens.analysis, bs)
            push!(fit_error, fit_effective_mass(subsample, corr, 1, fitting_range, binmethod = binmethod, nboot = nboot)[2])
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end

function plot_effective_mass(ens, corr, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = false)
    correlator_timeslices = eachindex(ens.analysis[1, corr])
    if plotting_range == :all
        plotting_range = (first(correlator_timeslices) + 1):last(correlator_timeslices)
    end
    μ, σ = bootstrap_effective_mass(ens.analysis, corr, binsize, binmethod=binmethod, nboot=nboot, range = plotting_range)
    plot_func = _bang ? plot! : plot
    plot_func(plotting_range, μ, yerr=σ)
end
function plot_effective_mass!(ens, corr, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = true)
    plot_effective_mass(ens, corr, binsize, plotting_range, binmethod = binmethod, nboot = nboot, _bang = _bang)
end

function plot_effective_mass_ratio(ens, numerator, denominator, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = false)
    μ, σ = bootstrap_effective_mass_ratio(ens.analysis, numerator, denominator, binsize, binmethod=binmethod, nboot=nboot)
    if plotting_range == :all
        plotting_range = eachindex(μ)
    end
    plot_func = _bang ? plot! : plot
    plot_func(plotting_range, μ[plotting_range], yerr=σ[plotting_range])
    xlabel!(L"$T$")
    ylabel!(String(numerator) * "/" * String(denominator) )
end
function plot_effective_mass_ratio!(ens, numerator, denominator, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = true)
    plot_effective_mass_ratio(ens, numerator, denominator, binsize, plotting_range, binmethod = binmethod, nboot = nboot, _bang = _bang)
end

function fit_effective_mass_ratio_naive(analysis::DataFrame, numerator, denominator, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = bootstrap_effective_mass_ratio(analysis, numerator, denominator, binsize, binmethod=binmethod, nboot=nboot)
    fit = fit_const(fitting_range, μ, σ)
    return [only(fit.param), only(stderror(fit))]
end

function fit_effective_mass_ratio(analysis::DataFrame, numerator, denominator, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    res = [@spawn try
                      subsample = get_subsample(analysis, binsize, method=binmethod)
                      μ, σ = bootstrap_effective_mass_ratio(subsample, numerator, denominator, 1, binmethod=binmethod, nboot=nboot)
                      only(fit_const(fitting_range, μ, σ).param)
           catch e
               missing
               end for i in 1:nboot]
    fit = fetch.(res)
    nfailures = sum(ismissing.(fit))
    if nfailures > 0
        @info "$nfailures bad bootstrap samples ($(100*nfailures/nboot)%)"
    end
    fit = [i for i in fit if !ismissing(i)]
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function tune_effective_mass_ratio_fit(ens, numerator, denominator, binsizes, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(ens.analysis, bs)
            push!(fit_error, fit_effective_mass_ratio(subsample, numerator, denominator, 1, fitting_range, binmethod = binmethod, nboot = nboot)[2])
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end
function tune_effective_mass_ratio_fit_naive(ens, numerator, denominator, binsizes, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(ens.analysis, bs)
            push!(fit_error, fit_effective_mass_ratio_naive(subsample, numerator, denominator, 1, fitting_range, binmethod = binmethod, nboot = nboot)[2])
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end

function plot_effective_mass_ratio_fit(ens, numerator, denominator, binsize, fitting_range, plotting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = fit_effective_mass_ratio(ens.analysis, numerator, denominator, binsize, fitting_range, binmethod=binmethod, nboot=nboot)
    @info μ, σ
    plot_effective_mass_ratio(ens, numerator, denominator, binsize, plotting_range, binmethod = binmethod, nboot = nboot)
    plot_const!(fitting_range, μ, σ)
    xlabel!(L"$T$")
    ylabel!(String(numerator) * "/" * String(denominator) )
end

function plot_effective_mass_ratio_fit_naive(ens, numerator, denominator, binsize, fitting_range, plotting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = fit_effective_mass_ratio_naive(ens.analysis, numerator, denominator, binsize, fitting_range, binmethod=binmethod, nboot=nboot)
    @info μ, σ
    plot_effective_mass_ratio(ens, numerator, denominator, binsize, plotting_range, binmethod = binmethod, nboot = nboot)
    plot_const!(fitting_range, μ, σ)
    xlabel!(L"$T$")
    ylabel!(String(numerator) * "/" * String(denominator) )
end

function plot_pcac_mass(ens, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = false)
    μ, σ = bootstrap_effective_pcac(ens.analysis, binsize, binmethod=binmethod, nboot=nboot)
    if plotting_range == :all
        plotting_range = eachindex(μ)
    end
    plot_func = _bang ? plot! : plot
    plot_func(plotting_range, μ[plotting_range], yerr=σ[plotting_range])
end
function plot_pcac_mass!(ens, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = true)
    plot_pcac_mass(ens, binsize, plotting_range, binmethod = binmethod, nboot = nboot, _bang = _bang)
end

function plot_fps(ens, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = false)
    L = last(ens.global_metadata[:geometry])
    T = first(ens.global_metadata[:geometry])
    μ, σ = bootstrap_effective_fps(ens.analysis, T, L, binsize, binmethod=binmethod, nboot=nboot)
    if plotting_range == :all
        plotting_range = eachindex(μ)
    end
    plot_func = _bang ? plot! : plot
    plot_func(plotting_range, μ[plotting_range], yerr=σ[plotting_range])
end
function plot_fps!(ens, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = true)
    plot_fps(ens, binsize, plotting_range, binmethod = binmethod, nboot = nboot, _bang = _bang)
end

function plot_gps(ens, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = false)
    L = last(ens.global_metadata[:geometry])
    T = first(ens.global_metadata[:geometry])
    μ, σ = bootstrap_effective_gps(ens.analysis, L, binsize, binmethod=binmethod, nboot=nboot)
    if plotting_range == :all
        plotting_range = eachindex(μ)
    end
    plot_func = _bang ? plot! : plot
    plot_func(plotting_range, μ[plotting_range], yerr=σ[plotting_range])
end
function plot_gps!(ens, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = true)
    plot_gps(ens, binsize, plotting_range, binmethod = binmethod, nboot = nboot, _bang = _bang)
end


function plot_effective_mass_fit(ens, corr, binsize, fitting_range, plotting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = fit_effective_mass(ens, corr, binsize, fitting_range, binmethod=binmethod, nboot=nboot)
    @info μ, σ
    plot_effective_mass(ens, corr, binsize, plotting_range, binmethod = binmethod, nboot = nboot)
    plot_const!(fitting_range, μ, σ)
    xlabel!(L"$T$")
    ylabel!(String(corr))
end
    
function plot_pcac_fit(ens, binsize, fitting_range, plotting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = fit_pcac_mass(ens, binsize, fitting_range, binmethod=binmethod, nboot=nboot)
    @info μ, σ
    plot_pcac_mass(ens, binsize, plotting_range, binmethod = binmethod, nboot = nboot)
    plot_const!(fitting_range, μ, σ)
end

function plot_gps_fit(ens, binsize, fitting_range, plotting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = fit_gps(ens, binsize, fitting_range, binmethod=binmethod, nboot=nboot)
    @info μ, σ
    plot_gps(ens, binsize, plotting_range, binmethod = binmethod, nboot = nboot)
    plot_const!(fitting_range, μ, σ)
end

function plot_fps_fit(ens, binsize, fitting_range, plotting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    μ, σ = fit_fps(ens, binsize, fitting_range, binmethod=binmethod, nboot=nboot)
    @info μ, σ
    plot_fps(ens, binsize, plotting_range, binmethod = binmethod, nboot = nboot)
    plot_const!(fitting_range, μ, σ)
end

function plot_dS(ens)
    plot(ens.analysis.dS, title = L"$dS$ History", xlabel = "Configuration")
end

function plot_maxeig(ens)
    range = ens.analysis[:, :confno]
    if length(Set(ens.analysis.maxeig)) == 1
        @error "No MaxEig data."
        return nothing
    else
        plot(range, ens.analysis.maxeig, title = "Maximum Eigenvalue", xlabel = "Configuration")
    end
end

function plot_mineig(ens)
    range = ens.analysis[:, :confno]
    if length(Set(ens.analysis.mineig)) == 1
        @error "No MinEig data."
        return nothing
    else
        plot(range, ens.analysis.mineig, title = "Minimum Eigenvalue", xlabel = "Configuration")
    end
end

function plot_adjoint_polyakov_hist(ens)
    data = hcat(ens.analysis.adjoint_polyakov...)
    histogram([data[i, :] for i in 1:4], layout = 4, plot_title = "Adjoint Polyakov", label = ["0" "1" "2" "3"])
end

function plot_adjoint_polyakov(ens)
    range = ens.analysis[:, :confno]
    data = hcat(ens.analysis.adjoint_polyakov...)
    plot(range, [data[i, :] for i in 1:4], layout = 4, plot_title = "Adjoint Polyakov", label = ["0" "1" "2" "3"])
end

function plot_fundamental_polyakov_hist(ens)
    data = real.(hcat(ens.analysis.fundamental_polyakov...))
    histogram([data[i, :] for i in 1:4], layout = 4, plot_title = "Fundamental Polyakov (Real Part)", label = ["0" "1" "2" "3"])
end

function plot_fundamental_polyakov(ens)
    range = ens.analysis[:, :confno]
    data = real.(hcat(ens.analysis.fundamental_polyakov...))
    plot(range, [data[i, :] for i in 1:4], layout = 4, plot_title = "Fundamental Polyakov (Real Part)", label = ["0" "1" "2" "3"])
end
