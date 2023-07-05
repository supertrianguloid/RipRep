using Plots
using Statistics
using Roots
using LsqFit
using DataFrames

include("parser.jl")
include("utilities.jl")

TITLE_FONT_SIZE = 13
NBOOT_DEFAULT = 50

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
    pcac = [OffsetVector{Float64, Vector{Float64}}[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        data = get_subsample(df, binsize, method=binmethod)
        push!(pcac[Threads.threadid()], effective_pcac(data))
    end
    pcac = reduce(vcat, pcac)
    return [mean(pcac), std(pcac)]
end

function bootstrap_effective_pcac_distribution(df::DataFrame, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    pcac = [OffsetVector{Float64, Vector{Float64}}[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        data = get_subsample(df, binsize, method=binmethod)
        push!(pcac[Threads.threadid()], effective_pcac(data))
    end
    return pcac
end

function fit_pcac_mass(ens, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    fit = [Float64[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(ens.analysis, binsize, method=binmethod)
        μ, σ = bootstrap_effective_pcac(subsample, 1, binmethod=binmethod, nboot=nboot)
        push!(fit[Threads.threadid()], only(fit_const(fitting_range, μ, σ).param))
    end
    fit = reduce(vcat, fit)

    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function effective_mass(correlator, T)
    meffs = Float64[]
    eq(m, τ) = h(T, τ - 1, 0, m)/h(T, τ, 0, m) - correlator[τ - 1]/correlator[τ]
    for τ in eachindex(correlator)[1:end]
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
    gps = [Vector{Float64}[] for _ in 1:Threads.nthreads()]
    T = length(df[1, :g5])
    
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        c = mean(subsample[:, :g5_folded])
        push!(gps[Threads.threadid()], effective_gps(c, T, L))
    end
    gps = reduce(vcat, gps)
    return [mean(gps), std(gps)]
end

function fit_gps(ens, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    fit = [Float64[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(ens.analysis, binsize, method=binmethod)
        μ, σ = bootstrap_effective_gps(subsample, L, 1, binmethod=binmethod, nboot=nboot)
        push!(fit[Threads.threadid()], only(fit_const(fitting_range, μ, σ).param))
    end
    fit = reduce(vcat, fit)
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function effective_fps(df::DataFrame, T, L)
    meff = effective_mass(mean(df[:, :g5_folded]), T)[1:end]
    gps = effective_gps(mean(df[:, :g5_folded]), T, L)[1:end]
    pcac = effective_pcac(df)
    return 2 .* (pcac./(meff.^2)).*gps
end

function bootstrap_effective_fps(df::DataFrame, T, L, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    fps = [OffsetArray{Float64}[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        push!(fps[Threads.threadid()], effective_fps(subsample, T, L))
    end
    fps = reduce(vcat, fps)
    return [mean(fps), std(fps)]
end

function fit_fps(ens, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    L = last(ens.global_metadata[:geometry])
    T = first(ens.global_metadata[:geometry])
    fit = [Float64[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(ens.analysis, binsize, method=binmethod)
        μ, σ = bootstrap_effective_fps(subsample, T, L, 1, binmethod=binmethod, nboot=nboot)
        push!(fit[Threads.threadid()], only(fit_const(fitting_range, μ, σ).param))
    end
    fit = reduce(vcat, fit)
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function bootstrap_effective_mass(df::DataFrame, corr, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT, range=:all)
    meff_boot = [Vector{Vector{Float64}}() for _ in 1:Threads.nthreads()]
    
    T = length(df[1, :g5])
    
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        c = mean(subsample[:, corr])
        if range != :all
            c = c[range]
        end
        
        effective_masses = effective_mass(c, T)
        push!(meff_boot[Threads.threadid()], effective_masses)
    end
    meff_boot = reduce(vcat, meff_boot)
    μ = mean(meff_boot)
    σ = std(meff_boot)
    return [μ, σ]
end

function bootstrap_effective_mass_ratio(df::DataFrame, corr_numerator, corr_denominator, binsize; binmethod = :randomsample, nboot = NBOOT_DEFAULT, range=:all)
    meff_boot = [Vector{Float64}[] for _ in 1:Threads.nthreads()]
    
    T = length(df[1, :g5])
    
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(df, binsize, method=binmethod)
        numerator = mean(subsample[:, corr_numerator])
        denominator = mean(subsample[:, corr_denominator])
        if range != :all
            numerator = numerator[range]
            denominator = denominator[range]
        end
        
        effective_masses = effective_mass_ratio(numerator, denominator, T)
        push!(meff_boot[Threads.threadid()], effective_masses)
    end
    meff_boot = reduce(vcat, meff_boot)
    μ = mean(meff_boot)
    σ = std(meff_boot)
    return [μ, σ]
end

function fit_effective_mass(ens, corr, binsize, fitting_range; binmethod = :randomsample, nboot = NBOOT_DEFAULT)
    fit = [Float64[] for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:nboot
        subsample = get_subsample(ens.analysis, binsize, method=binmethod)
        μ, σ = bootstrap_effective_mass(subsample, corr, 1, binmethod=binmethod, nboot=nboot)
        push!(fit[Threads.threadid()], only(fit_const(fitting_range, μ, σ).param))
    end
    fit = reduce(vcat, fit)
    μ = mean(fit)
    σ = std(fit)
    return [μ, σ]
end

function plot_effective_mass(ens, corr, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = false)
    μ, σ = bootstrap_effective_mass(ens.analysis, corr, binsize, binmethod=binmethod, nboot=nboot)
    if plotting_range == :all
        plotting_range = eachindex(μ)
    end
    plot_func = _bang ? plot! : plot
    plot_func(plotting_range, μ[plotting_range], yerr=σ[plotting_range])
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
end
function plot_effective_mass_ratio!(ens, numerator, denominator, binsize, plotting_range = :all; binmethod = :randomsample, nboot = NBOOT_DEFAULT, _bang = true)
    plot_effective_mass(ens, numerator, denominator, binsize, plotting_range, binmethod = binmethod, nboot = nboot, _bang = _bang)
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
