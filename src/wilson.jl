using Base: nothing_sentinel, full_va_len
using Plots
using LsqFit
using Base.Threads: @spawn
using LaTeXStrings

NBOOT_DEFAULT = 20

include("parser.jl")
include("utilities.jl")

function plot_w(wf, range = :all; binsize = 1, nboot = 1000, _bang = false, title = "", sym = true)
    plot_func = _bang ? plot! : plot
    data = wf.analysis
    if range == :all
        range = data[1, :t][2:end-1]
    end
    if(sym)
        data = data[:, :Wsym]
    else
        data = data[:, :W]
    end

    index = _wf_time_to_index_w(wf.analysis, range[1]):_wf_time_to_index_w(wf.analysis, range[end])
    plot_func(index.*wf.metadata[:dt], mean(data)[index], yerr = standard_error(data, binsize=binsize, nboot=nboot)[index], title=title, legend=false)
    xlabel!(L"t")
    ylabel!(L"W(t)")
end

function plot_w!(wf, range = :all, sym = true)
    plot_w(wf, range, _bang = true, sym = sym)
end

function plot_tc(wf; title="")
    a = findall(x -> x == 1, mean(wf.data[:,:W]) .> 1)
    if !isempty(a)
        t = first(a)
    else
        t = length(wf.data.t[1])÷2
    end
    plot(1:nrow(wf.data), [i[t] for i in wf.data[:, :TC]], title=title, label=latexstring("t = ", wf.data[1,:t][t]))
end

function plot_tc_hist(wf; title="")
    a = findall(x -> x == 1, mean(wf.data[:,:W]) .> 1)
    if !isempty(a)
        t = first(a)
    else
        t = length(wf.data.t[1])÷2
    end
    histogram([i[t] for i in wf.data[:, :TC]], title=title, label=latexstring("t = ", wf.data[1,:t][t]))
end

function plot_t2e(wf, range = :all; binsize = 1, nboot=1000, title="")
    data = wf.analysis
    if range == :all
        range = data[1, :t]
    end
    index = _wf_time_to_index(wf.analysis, range[1]):_wf_time_to_index(wf.analysis, range[end])
    plot((index .- 1).*wf.metadata[:dt], mean(data[:, :t2E])[index], yerr = standard_error(data[:, :t2E], binsize=binsize, nboot=nboot)[index], title=title, legend=false)
    xlabel!("\$t\$")
    ylabel!("\$t^2E(t)\$")
end

function _wf_time_to_index(analysis::DataFrame, time)
    return only(findall(analysis.t[1] .≈ time))
end

function _wf_time_to_index_w(analysis::DataFrame, time)
    return only(findall(analysis.t[1] .≈ time)) - 1
end

function _find_t0(analysis::DataFrame, data, window, dt; nboot = 100, ref = 1.0)
    n1 = _wf_time_to_index(analysis, window[1])
    n2 = _wf_time_to_index(analysis, window[2])
    res = [@spawn begin
                    s = rand(1:nrow(data), nrow(data))
                    σ² = var(data[s, :t2E])
                    mu = mean(data[s, :t2E])
                    @. model(x, params) = params[1]*x + params[2]
                    fit = curve_fit(model, (n1:n2).*dt, mu[n1:n2], (1 ./σ²)[n1:n2], [0.0, 0.0])
                    [fit.param[1], fit.param[2]]
                  end for i in 1:nboot]
    a = fetch.(res)
    m = vcat(a'...)[:,1]
    c = vcat(a'...)[:,2]
    
    t = reference_time(ref, mean(c), var(c), mean(m), var(m), cov(m, c))
    return t
end

function _find_w0(wf, window, dt; binsize = 1, nboot = 100, ref = 1.0, sym = true)
    n1 = _wf_time_to_index_w(wf, window[1])
    n2 = _wf_time_to_index_w(wf, window[2])
    res = [@spawn begin
                    sample = get_subsample(wf.analysis, binsize)
                    if(sym)
                        σ² = var(sample[:, :Wsym])
                        mu = mean(sample[:, :Wsym])
                    else
                        σ² = var(sample[:, :W])
                        mu = mean(sample[:, :W])
                    end
                    @. model(x, params) = params[1]*x + params[2]
                    fit = curve_fit(model, (n1:n2).*dt, mu[n1:n2], (1 ./σ²)[n1:n2], [0.0, 0.0])
                    [fit.param[1], fit.param[2]]
                    end for i in 1:nboot]
    a = fetch.(res)
    m = vcat(a'...)[:,1]
    c = vcat(a'...)[:,2]
    w = reference_time(ref, mean(c), var(c), mean(m), var(m), cov(m, c))
    return [sqrt(w[1]), w[2]/(2*sqrt(w[1]))]
end

function _find_w0_fullbootstrap(analysis::DataFrame, window, dt; binsize = 1, nboot = 100, ref = 1.0, sym = true)
    n1 = _wf_time_to_index_w(analysis, window[1])
    n2 = _wf_time_to_index_w(analysis, window[2])
    res = [@spawn begin
                    sample = get_subsample(analysis, binsize)
                    if(sym)
                        σ² = var(sample[:, :Wsym])
                        mu = mean(sample[:, :Wsym])
                    else
                        σ² = var(sample[:, :W])
                        mu = mean(sample[:, :W])
                    end
                    @. model(x, params) = params[1]*x + params[2]
                    fit = curve_fit(model, (n1:n2).*dt, mu[n1:n2], (1 ./σ²)[n1:n2], [0.0, 0.0])
                    m, c = fit.param
                    sqrt((ref - c)/m)
                end for i in 1:nboot]
    w = fetch.(res)
    return [mean(w), std(w)]
end

function find_w0(analysis::DataFrame, dt, window; binsize = 1, nboot = 100, ref = 1.0, sym = true, fullbootstrap = true)
    if(fullbootstrap)
        return _find_w0_fullbootstrap(analysis, window, dt, binsize = binsize, nboot = nboot, ref = ref, sym = sym)
    else
        return _find_w0(analysis, window, dt, binsize = binsize, nboot = nboot, ref = ref, sym = sym)
    end
end

function find_w0(wf::WilsonFlow, window; binsize = 1, nboot = 100, ref = 1.0, sym = true, fullbootstrap = true)
    return find_w0(wf.analysis, wf.metadata[:dt], window, binsize=binsize, nboot=nboot, ref=ref, sym=sym, fullbootstrap=fullbootstrap)
end

function find_t0(wf, window; nboot = 100, ref = 1.0)
    return _find_t0(wf, wf.analysis, window, wf.metadata[:dt], nboot = nboot, ref = ref)
end

function reference_time(yref, c, cvar, m, mvar, cov = 0)
    return (yref - c)/m, sqrt(
        (1/m)^2 * cvar + 
        (((c - yref)^2)/m^4) * mvar + 
        2*cov*(yref - c)/(m^3)
    )
end

function auto_w0(analysis::DataFrame, dt; binsize = 1, eitherside = 2, nboot = 100, ref = 1, sym=true, fullbootstrap=true)
    if(sym)
        data = analysis.Wsym
    else
        data = analysis.W
    end
    rightpoint = findfirst(>(ref), mean(data))
    window = 0
    try
        fitrange = (rightpoint-eitherside):((rightpoint - 1) + eitherside)
        a = mean(data)[fitrange]
        window = analysis.t[1][fitrange]
    catch e
        @error "W never reaches the reference value." ref
        return nothing
    end
    return find_w0(analysis, dt, window, binsize = binsize, nboot = nboot, ref = ref, sym = sym, fullbootstrap = fullbootstrap)
end

function auto_w0(wf::WilsonFlow; binsize = 1, eitherside = 2, nboot = NBOOT_DEFAULT, ref = 1, sym=true, fullbootstrap=true)
    return auto_w0(wf.analysis, wf.metadata[:dt], binsize = binsize, eitherside = eitherside, nboot = nboot, ref = ref, sym=sym, fullbootstrap=fullbootstrap)
end
    

function tune_binsize_auto_w0(wf::WilsonFlow, binsizes; eitherside = 2, ref = 1, sym = true, nboot = NBOOT_DEFAULT)
    y = []
    yerr = []
    for bs in binsizes
        fit_error = []
        for i in 1:nboot
            subsample = get_subsample(wf.analysis, bs)
            try
                e = auto_w0(subsample, wf.metadata[:dt], binsize = 1, eitherside = eitherside, nboot = nboot, ref = ref, sym = sym)[2]
                push!(fit_error, e)
            catch e
                
            end
        end
        push!(y, mean(fit_error))
        push!(yerr, std(fit_error))
    end
    plot(binsizes, y, yerr=yerr, xticks = binsizes)
end
