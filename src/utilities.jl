function _ensemble_to_latex_string(ens)
    return "\$\\beta = " * string(only(ens.global_metadata.β)) * ",\\ C_{SW} = " * string(only(ens.global_metadata.Cˢʷ)) * ",\\ m = " * string(only(ens.global_metadata.Mass)) * ",\\ V = " * string(only(ens.global_metadata.T)) * "\\times" * string(only(ens.global_metadata.L)) * "^3,\\ \\mathrm{" * only(ens.global_metadata.SimulationType) * "}\$"
end

function _only_match(regex, output)
    m = match.(regex, output)
    return only(only(m[m .!= nothing]).captures)
end