using Plots
using LsqFit

include("parser.jl")
include("utilities.jl")

function plot_w(wf, range = :all; _bang = false, title = "")
    plot_func = _bang ? plot! : plot
    data = wf.analysis
    if range == :all
        range = data[1, :t][2:end-1]
    end
    index = _wf_time_to_index_w(wf, range[1]):_wf_time_to_index_w(wf, range[end])
    plot_func(index.*wf.metadata[:dt], mean(data[:, :W])[index], yerr = _bootstrap(data[:, :W])[index], title=title, legend=false)
    xlabel!("\$t\$")
    ylabel!("\$W(t)\$")
end

function plot_w!(wf, range = :all)
    plot_w(wf, range, _bang = true)
end

function plot_tc(wf; title="")
    a = findall(x -> x == 1, mean(wf.data[:,:W]) .> 1)
    if !isempty(a)
        t = first(a)
    else
        t = length(wf.data.t[1])÷2
    end
    plot(1:nrow(wf.data), [i[t] for i in wf.data[:, :TC]], title=title, label="t = "*string(wf.data[1,:t][t]))
end

function plot_tc_hist(wf; title="")
    a = findall(x -> x == 1, mean(wf.data[:,:W]) .> 1)
    if !isempty(a)
        t = first(a)
    else
        t = length(wf.data.t[1])÷2
    end
    histogram([i[t] for i in wf.data[:, :TC]], title=title, label="t = "*string(wf.data[1,:t][t]))
end

function plot_t2e(wf, range = :all; nboot = 1000, title="")
    data = wf.analysis
    if range == :all
        range = data[1, :t]
    end
    index = _wf_time_to_index(wf, range[1]):_wf_time_to_index(wf, range[end])
    plot(index.*wf.metadata[:dt], mean(data[:, :t2E])[index], yerr = _bootstrap(data[:, :t2E], nboot)[index], title=title, legend=false)
    xlabel!("\$t\$")
    ylabel!("\$t^2E(t)\$")
end

function _bootstrap(data, nboot = 1000)
    l = length(data)
    bs = []
    for i ∈ 1:nboot
        push!(bs, mean(data[rand(1:l, l)]))
    end
    return std(bs)
end

function _wf_time_to_index(wf::WilsonFlow, time)
    return only(findall(wf.analysis.t[1] .≈ time))
end

function _wf_time_to_index_w(wf::WilsonFlow, time)
    return only(findall(wf.analysis.t[1] .≈ time)) - 1
end

function _find_t0(wf, data, window, dt; nboot = 100, ref = 1.0)
    n1 = _wf_time_to_index(wf, window[1])
    n2 = _wf_time_to_index(wf, window[2])
    m = []
    c = []
    for i in 1:nboot
        s = rand(1:nrow(data), nrow(data))
        σ² = var(data[s, :t²E])
        mu = mean(data[s, :t²E])
        @. model(x, params) = params[1]*x + params[2]
        fit = curve_fit(model, (n1:n2).*dt, mu[n1:n2], (1 ./σ²)[n1:n2], [0.0, 0.0])
        push!(m, fit.param[1])
        push!(c, fit.param[2])
    end
    t = reference_time(ref, mean(c), var(c), mean(m), var(m), cov(m, c))
    return t
end

function _find_w0(wf, data, window, dt; nboot = 100, ref = 1.0)
    n1 = _wf_time_to_index_w(wf, window[1])
    n2 = _wf_time_to_index_w(wf, window[2])
    m = []
    c = []
    for i in 1:nboot
        s = rand(1:nrow(data), nrow(data))
        σ² = var(data[s, :W])
        mu = mean(data[s, :W])
        @. model(x, params) = params[1]*x + params[2]
        fit = curve_fit(model, (n1:n2).*dt, mu[n1:n2], (1 ./σ²)[n1:n2], [0.0, 0.0])
        push!(m, fit.param[1])
        push!(c, fit.param[2])
    end
    w = reference_time(ref, mean(c), var(c), mean(m), var(m), cov(m, c))
    return sqrt(w[1]), w[2]/(2*sqrt(w[1]))
end

function find_w0(wf, window; nboot = 100, ref = 1.0)
    return _find_w0(wf, wf.analysis, window, wf.metadata.dt, nboot = nboot, ref = ref)
end

function find_t0(wf, window; nboot = 100, ref = 1.0)
    return _find_t0(wf, wf.analysis, window, wf.metadata.dt, nboot = nboot, ref = ref)
end

function error_on_error_w0(wf, window; nboot = 100, ref = 1.0)
    w0 = []
    for i in 1:nboot
        push!(w0, _find_w0(wf, wf.analysis[rand(1:nrow(wf.analysis), nrow(wf.analysis)), :], window, wf.metadata.dt, nboot=nboot, ref=ref))
    end
    w0 = collect.(w0)
    return hcat(mean(w0)..., std(w0)...)
end

function error_on_error_t0(wf, window; nboot = 100, ref = 1.0)
    t0 = []
    for i in 1:nboot
        push!(t0, _find_t0(wf, wf.analysis[rand(1:nrow(wf.analysis), nrow(wf.analysis)), :], window, wf.metadata.dt, nboot=nboot, ref=ref))
    end
    t0 = collect.(t0)
    return hcat(mean(t0)..., std(t0)...)
end

function wf_plot_error_on_error(wf, window, binrange, type; nboot = 100, ref = 1.0, method = :equal)
    if type == :w0
        l1, l2 = mean(wf.data[:, :W])[[_wf_time_to_index_w(wf, window[1]), _wf_time_to_index_w(wf, window[2])]]
    elseif type == :t0
        l1, l2 = mean(wf.data[:, :t²E])[[_wf_time_to_index(wf, window[1]), _wf_time_to_index(wf, window[2])]]
    end
    if l1 > ref || l2 < ref
        @error "Wrong bounds" * string(window)
        return
    end
    res = []
    for i in binrange
        thermalise_bin!(wf, 1, i, method = method)
        if type == :w0
            push!(res, error_on_error_w0(wf, window, nboot = nboot, ref = ref))
        elseif type == :t0
            push!(res, error_on_error_t0(wf, window, nboot = nboot, ref = ref))
        end
    end
    res = vcat(res...)
    plot(binrange, res[:,3], yerr = res[:,4])
end

function reference_time(yref, c, cvar, m, mvar, cov = 0)
    return (yref - c)/m, sqrt(
        (1/m)^2 * cvar + 
        (((c - yref)^2)/m^4) * mvar + 
        2*cov*(yref - c)/(m^3)
    )
end

#function find_w0(wf, window, ref = 1.0)
#end
