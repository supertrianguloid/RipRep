using DrWatson
using Plots
using Statistics
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

function _plot_plaquette(run_no, data)
    plot(1:nrow(data), data[:, :plaquette], label="Plaquette", title = GLOBAL_SIMS[run_no, :path])
    xlabel!("Conf #")
end

function plot_correlator(correlator)
    plot(0:length(mean(correlator)) - 1, parent(mean(correlator)), yerr = parent(std(correlator, corrected=true)/sqrt(length(correlator))))
    xlabel!("Imaginary Time")
end

function fit_cosh(data, correlator, trange, thermalisation, binsize, nstates, bs, mode)
    
end