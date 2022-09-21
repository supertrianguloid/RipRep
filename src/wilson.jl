using Plots

theme(:dracula)
default(size = (900, 700))

include("parser.jl")
include("utilities.jl")

@info "RipRep Wilson flow code activated."

function plot_t2e(wf)
    data = wf.analysis
    plot(data[1, :t], mean(data[:, :t²E]), yerr = _bootstrap(data[:, :t²E]))
    xlabel!("\$t\$")
    ylabel!("\$t^2E\$")
end

function bin!(wf::WilsonFlow, binsize, method = :equal)
    if binsize == 1
        return
    end
    data = wf.data
    offset = nrow(data) % binsize
    analysis = data[offset:end, :]
    newlen = nrow(analysis) ÷ binsize
    if method == :equal
        wf.analysis = analysis[[1 + (i - 1)*binsize for i in 1:newlen], :]
        return
    end
end

function _bootstrap(data, n = 1000)
    l = length(data)
    bs = []
    for i ∈ 1:n
        push!(bs, mean(data[rand(1:l, l)]))
    end
    return std(bs)
end