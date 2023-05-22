using Statistics

function thermalise_bin!(data_struct, thermsize, binsize; method = :equal)
    dtype = typeof(data_struct)
    data = data_struct.data[thermsize:end,:]
    if binsize == 1
        data_struct.analysis = data
        return
    end

    offset = nrow(data) % binsize + 1
    data = data[offset:end, :]
    newlen = nrow(data) ÷ binsize

    if method == :equal
        data_struct.analysis = data[[1 + (i - 1)*binsize for i in 1:newlen], :]
        return
    elseif method == :average
        if dtype == Ensemble
            sample = select(data, Not([:runno, :time, :completed, :accepted]))
        elseif dtype == WilsonFlow
            sample = select(data, Not([:confno]))
        end
        df = DataFrame()
        for i in 1:newlen
            index = ((i-1)*binsize + 1):(i*binsize)
            row = [[x] for x in mean.(eachcol(sample[index,:]))]
            push!(df, DataFrame(row, names(sample))[1,:])
        end
        data_struct.analysis = df
    else
        @error "Unknown bin method" method
    end

end

function propagate_product(x, y)
    return [x[1]*y[1], sqrt((x[2]/x[1])^2 + (y[2]/y[1])^2)]
end

function propagate_square(x)
    return [x[1]^2, (2*x[1]*x[2])]
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

function _only_match(regex, output)
    m = match.(regex, output)
    return only(only(m[m .!= nothing]).captures)
end
function _ensemble_to_latex_string(ens)
    return "\$\\beta = " * string(only(ens.global_metadata[:β])) * ",\\ C_{SW} = " * string(only(ens.global_metadata[:csw])) * ",\\ m = " * string(only(ens.global_metadata[:m0])) * ",\\ V = " * string(only(ens.global_metadata[:geometry][0])) * "\\times" * string(only(ens.global_metadata[:geometry][1])) * "\\times" * string(only(ens.global_metadata[:geometry][2])) * "\\times" * string(only(ens.global_metadata[:geometry][3])) * "\$"
end
