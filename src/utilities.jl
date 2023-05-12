function _ensemble_to_latex_string(ens)
    return "\$\\beta = " * string(only(ens.global_metadata.β)) * ",\\ C_{SW} = " * string(only(ens.global_metadata.Cˢʷ)) * ",\\ m = " * string(only(ens.global_metadata.Mass)) * ",\\ V = " * string(only(ens.global_metadata.T)) * "\\times" * string(only(ens.global_metadata.L)) * "^3\$"
end

function _only_match(regex, output)
    m = match.(regex, output)
    return only(only(m[m .!= nothing]).captures)
end

function d(a; h = 1)
    start = a.offsets[1] + 1
    ed = eachindex(a)[end]
    deriv = []
    for t in Iterators.take(eachindex(a), length(a)-2)
        push!(deriv, a[t+2] - a[t])
    end
    return OffsetArray(deriv, (start+1):ed-1)./(2*h)
end

function propagate_product(x, y)
    return [x[1]*y[1], sqrt((x[2]/x[1])^2 + (y[2]/y[1])^2)]
end

function propagate_square(x)
    return [x[1]^2, (2*x[1]*x[2])]
end
