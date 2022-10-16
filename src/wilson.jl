using Plots
using LsqFit

theme(:dracula)
default(size = (900, 700))

include("parser.jl")
include("utilities.jl")

@info "RipRep Wilson flow code activated."

function plot_w(wf, range = :all)
    data = wf.analysis
    if range == :all
        range = data[1, :t][2:end-1]
    end
    index = _wf_time_to_index_w(wf, range[1]):_wf_time_to_index_w(wf, range[end])
    plot(index.*wf.metadata.dt, mean(data[:, :W])[index], yerr = _bootstrap(data[:, :W])[index])
    xlabel!("\$t\$")
    ylabel!("\$W(t)\$")
end

function plot_t2e(wf, range = :all, nboot = 1000)
    data = wf.analysis
    if range == :all
        range = data[1, :t]
    end
    index = _wf_time_to_index(wf, range[1]):_wf_time_to_index(wf, range[end])
    plot(index.*wf.metadata.dt, mean(data[:, :t²E])[index], yerr = _bootstrap(data[:, :t²E], nboot)[index])
    xlabel!("\$t\$")
    ylabel!("\$t^2E(t)\$")
end

function plot_t2e!(wf, range = :all)
    if range == :all
        range = data[1, :t]
    end
    data = wf.analysis
    plot!(range, mean(data[:, :t²E]), yerr = _bootstrap(data[:, :t²E]))
    xlabel!("\$t\$")
    ylabel!("\$t^2E\$")
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
    return only(findall(wf.analysis.t[1] .== time))
end

function _wf_time_to_index_w(wf::WilsonFlow, time)
    return only(findall(wf.analysis.t[1] .== time)) - 1
end

function find_t0(wf, window; bs = 100, ref = 1.0)
    n1 = _wf_time_to_index(wf, window[1])
    n2 = _wf_time_to_index(wf, window[2])
    m = []
    c = []
    for i in 1:bs
        s = rand(1:nrow(wf.analysis), nrow(wf.analysis))
        sigma = std(wf.analysis[s, :t²E])
        mu = mean(wf.analysis[s, :t²E])
        @. model(x, params) = params[1]*x + params[2]
        f = curve_fit(model, (n1:n2).*wf.metadata.dt, mu[n1:n2], (1 ./sigma)[n1:n2], [0.0, 0.0]).param
        push!(m, f[1])
        push!(c, f[2])
    end
    return reference_time(ref, mean(c), std(c), mean(m), std(m))
end

function find_w0(wf, window; bs = 100, ref = 1.0)
    n1 = _wf_time_to_index_w(wf, window[1])
    n2 = _wf_time_to_index_w(wf, window[2])
    m = []
    c = []
    for i in 1:bs
        s = rand(1:nrow(wf.analysis), nrow(wf.analysis))
        sigma = std(wf.analysis[s, :W])
        mu = mean(wf.analysis[s, :W])
        @. model(x, params) = params[1]*x + params[2]
        f = curve_fit(model, (n1:n2).*wf.metadata.dt, mu[n1:n2], (1 ./sigma)[n1:n2], [0.0, 0.0]).param
        push!(m, f[1])
        push!(c, f[2])
    end
    w = reference_time(ref, mean(c), std(c), mean(m), std(m))
    return sqrt(w[1]), 1/(2*sqrt(w[2]))
end

function reference_time(yref, c, cerr, m, merr, cov = 0)
    return (yref - c)/m, sqrt(
        (1/m)^2 * (cerr^2) + 
        (((c - yref)^2)/m^4) * (merr^2) + 
        2*cov*(yref - c)/(m^3)
    )
end

#function find_w0(wf, window, ref = 1.0)
#end