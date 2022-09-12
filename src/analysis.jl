
function acceptance(ens)
    return sum(ens.data[:, :accepted])/nrow(ens.data)
end

function thermalise!(ens::Ensemble, ntherm)
    ens.analysis = ens.data[ntherm:end, :]
end

function thermalise(data::DataFrame, ntherm)
    return data[ntherm:end, :]
end

function bin(ens, binsize, method = :equal)
    data = ens.data
    if binsize == 1
        return
    end
    offset = nrow(data) % binsize
    data = thermalise(data, offset)
    newlen = nrow(data) ÷ binsize
    if method == :equal
        ens.analysis = data[[1 + (i - 1)*binsize for i in 1:newlen], :]
    end
end

function bin(data::DataFrame, binsize, method = :equal)
    if binsize == 1
        return
    end
    offset = nrow(data) % binsize
    data = thermalise(data, offset)
    newlen = nrow(data) ÷ binsize
    if method == :equal
        return data[[1 + (i - 1)*binsize for i in 1:newlen], :]
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

function plot_t2e(wf, binsize = 1)
    data = bin(wf.data, binsize)
    plot(data[1, :t], mean(data[:, :t²E]), yerr = _bootstrap(data[:, :t²E]))
    xlabel!("\$t\$")
    ylabel!("\$t^2E\$")
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
    plot(tmin:tmax - nstates*2, fits[:, 2], yerr = fits[:, 4], label = string(correlator) * " mass, \$\\tau_{max}\$ = $tmax", title =  _ensemble_to_latex_string(ens))
    xlabel!("Lower fitting range")
end