using DrWatson
@quickactivate "RipRep"
include(srcdir("parser.jl"))

function acceptance(data)
    return sum(data[:, :accepted])/nrow(data)
end

function plot_plaquette(data)

end

function fit_cosh(data, correlator, trange, thermalisation, binsize, nstates, bs, mode)
    
end