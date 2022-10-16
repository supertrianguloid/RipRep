
function acceptance(ens)
    return sum(ens.data[:, :accepted])/nrow(ens.data)
end

function thermalise(data::DataFrame, ntherm)
    return data[ntherm:end, :]
end

function plot_t2e(wf, binsize = 1)
    data = bin(wf.data, binsize)
    plot(data[1, :t], mean(data[:, :t²E]), yerr = _bootstrap(data[:, :t²E]))
    xlabel!("\$t\$")
    ylabel!("\$t^2E\$")
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