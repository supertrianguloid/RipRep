using Statistics

function thermalise!(data_struct, thermsize; method = :equal)
    data_struct.analysis = data_struct.data[thermsize:end,:]
end


# TODO: Fix this to stop offsetting the beginning. Don't need this anymore

function get_subsample(df::DataFrame, binsize; method = :randomsample)
    offset = nrow(df) % binsize + 1
    subsample = df[offset:end, :]
    newlen = nrow(subsample) ÷ binsize

    if method == :randomsample
        return subsample[rand(1:nrow(subsample), newlen), :]
    end
        

    # if method == :equal
    #     data_struct.analysis = data[[1 + (i - 1)*binsize for i in 1:newlen], :]
    #     return
    # elseif method == :average
    #     if dtype == Ensemble
    #         sample = select(data, Not([:runno, :time, :completed, :accepted]))
    #     elseif dtype == WilsonFlow
    #         sample = select(data, Not([:confno]))
    #     end
    #     df = DataFrame()
    #     for i in 1:newlen
    #         index = ((i-1)*binsize + 1):(i*binsize)
    #         row = [[x] for x in mean.(eachcol(sample[index,:]))]
    #         push!(df, DataFrame(row, names(sample))[1,:])
    #     end
    #     data_struct.analysis = df
    # else
    #     @error "Unknown bin method" method
    # end
end

function get_subsample(vec::Vector, binsize; method = :randomsample)
    offset = length(vec) % binsize + 1
    subsample = vec[offset:end]
    newlen = length(subsample) ÷ binsize

    if method == :randomsample
        return subsample[rand(1:length(subsample), newlen)]
    end
end

function standard_error(data; binsize = 1, nboot = 1000)
    bootstrap = []
    for i in 1:nboot
        sample = get_subsample(data, binsize)
        push!(bootstrap, mean(sample))
    end
    return std(bootstrap)
end

function propagate_product(x, y)
    return [x[1]*y[1], sqrt((x[2]/x[1])^2 + (y[2]/y[1])^2)]
end

function propagate_square(x)
    return [x[1]^2, (2*x[1]*x[2])]
end

function d(a; h = 1, last_point_antisymmetric = false)
    array_indices = eachindex(a)
    deriv = []
    for t in collect(array_indices)[2:end-1]
        push!(deriv, a[t+1] - a[t-1])
    end
    if last_point_antisymmetric
        push!(deriv, -2*a[end-1])
        return OffsetArray(deriv, (first(array_indices)+1):last(array_indices))./(2*h)
    else
        return OffsetArray(deriv, (first(array_indices)+1):last(array_indices)-1)./(2*h)
    end
end

function _corr_derivative(correlator; h = 1, last_point_antisymmetric = false)
    corr = []
    for c in correlator
        if c === missing
            push!(corr, missing)
            continue
        end
        push!(corr, d(c, h=h, last_point_antisymmetric=last_point_antisymmetric))
    end
    return corr
end

function fit_const(range, μ, σ)
    μ = μ[range]
    σ = σ[range]
    w = 1 ./ (σ.^2)
    return curve_fit((x, params) -> only(params) .+ 0*x, range, μ, w, [0.0])
end

function fit_line(range, μ, σ)
    μ = μ[range]
    σ = σ[range]
    # Fits y = c₁x + c₂
    w = 1 ./ (σ.^2)
    return curve_fit((x, params) -> params[1].*x .+ params[2], range, μ, w, [0.0, 0.0])
end

function plot_line(range, params; _bang = false)
    plot_func = _bang ? plot! : plot
    y = collect(range).*params[1] .+ params[2]
    plot_func(range, y)
end

function plot_line!(range, params; _bang = true)
    plot_line(range, params, _bang = true)
end

function plot_const(range, c, cerr = nothing; _bang = false)
    plot_func = _bang ? plot! : plot
    if cerr != nothing
        plot_func(range, repeat([c], length(range)), ribbon=(cerr, cerr))
    else
        plot_func(range, repeat([c], length(range)))
    end
end

function plot_const!(range, c, cerr = nothing)
    plot_const(range, c, cerr, _bang = true)
end

function plot_in_range(range, μ, σ; _bang = false)
    μ = μ[range]
    σ = σ[range]
    plot_func = _bang ? plot! : plot
    plot_func(range, μ, yerr = σ)
end

function plot_in_range!(range, μ, σ)
    plot_in_range(range, μ, σ, _bang=true)
end


function _only_match(regex, output)
    m = match.(regex, output)
    return only(only(m[m .!= nothing]).captures)
end

function _ensemble_to_latex_string(ens)
    return "\$\\beta = " * string(only(ens.global_metadata[:β])) * ",\\ C_{SW} = " * string(only(ens.global_metadata[:csw])) * ",\\ m = " * string(only(ens.global_metadata[:m0])) * ",\\ V = " * string(only(ens.global_metadata[:geometry][0])) * "\\times" * string(only(ens.global_metadata[:geometry][1])) * "\\times" * string(only(ens.global_metadata[:geometry][2])) * "\\times" * string(only(ens.global_metadata[:geometry][3])) * "\$"
end

function _fold(correlator; symm = true)
    if symm == true
        sign = 1
    else
        sign = -1
    end
    corr = []
    tmax = length(correlator[1])
    for c in correlator
        if c === missing
            push!(corr, missing)
            continue
        end
        v = [c[0]]
        for i in 1:((tmax÷2)-1)
            push!(v, 0.5*(c[i] + sign*c[tmax - i]))
        end
        push!(v, c[tmax÷2])
        push!(corr, OffsetArray(v, 0:tmax÷2))
    end
    return corr
end

function string_to_unitrange(str::String)
    arr = parse.(Int, split(str, ":"))
    return first(arr):last(arr)
end

function ensure_directory_exists(path)
    try
        mkpath(path)
    catch e
    end
end
