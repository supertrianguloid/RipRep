using Plots
using LsqFit
using DataFrames
using LaTeXStrings
using PyCall
using SciPy
include("utilities.jl")

function plot_data(data, x, y; β = :all, csw = :all, _bang = false, markersize = 1, markershape = :circle, ignore = [])
    plot_func = _bang ? scatter! : scatter
    x_vals, x_err, y_vals, y_err = get_data(data, x, y, β = β, csw = csw, ignore = ignore)

    plot_func(x_vals, y_vals, xlabel = String(x), ylabel = String(y), yerr = y_err, xerr = x_err, markersize=markersize, markershape = markershape, legend = false, label = L"$\beta = %$β$" )
end

function plot_data!(data, x, y; β = :all, csw = :all, _bang = true, markersize = 1, markershape = :circle, ignore = [])
    plot_data(data, x, y, β = β, csw = csw, _bang = _bang, markersize = markersize, markershape = markershape, ignore = ignore)
end

function get_data(data, x, y; β = :all, csw = :all, ignore = [], xlims = [-100000, 100000], ylims = [-100000, 100000])
    data = data[[i for i=1:nrow(data) if !(i in ignore)], :]
    filtered = DataFrames.filter([:csw, :β, x, y] => (csw_v, β_v, x_v, y_v) -> (csw_v == csw || csw == :all) && (β_v == β || β == :all) && !ismissing(x_v) && !isnothing(x_v) && !ismissing(y_v) && !isnothing(y_v) && y_v[1] > ylims[1] && y_v[1] < ylims[2] && x_v[1] > xlims[1] && x_v[1] < xlims[2], data)
    y_vals = filtered[:, y]
    y_err = nothing
    x_vals = filtered[:, x]
    x_err = nothing
    if isa(first(y_vals), Vector) && length(first(y_vals)) == 2
        y_err = [i[2] for i in y_vals ]
        y_vals = [i[1] for i in y_vals]
    end
    if isa(first(x_vals), Vector) && length(first(x_vals)) == 2
        x_err = [i[2] for i in x_vals]
        x_vals = [i[1] for i in x_vals]
    end
    return x_vals, x_err, y_vals, y_err
end

function linear_fit(data, x, y;  β = :all, csw = :all, _bang = false, markersize = 1, markershape = :circle, ignore = [], plotrange = :default, xlims = [-100000, 100000], ylims = [-100000, 100000])
    x_vals, x_err, y_vals, y_err = get_data(data, x, y, β = β, csw = csw, ignore = ignore, xlims = xlims, ylims = ylims)
    plot_data(data, x, y, β = β, csw = csw, _bang = _bang, markersize = markersize, markershape = markershape, ignore = ignore)
    range = [min(x_vals...), max(x_vals...)]
    cov = nothing
    params = nothing
    chisq = nothing
    if isnothing(x_err)
        fit = nothing
        if !isnothing(y_err)
            w = 1 ./ (y_err.^2)
            fit = curve_fit((x, params) -> params[1].*x .+ params[2], x_vals, y_vals, w, [0.0, 0.0])
        else
            fit = curve_fit((x, params) -> params[1].*x .+ params[2], x_vals, y_vals, [0.0, 0.0])
        end
        
        cov = estimate_covar(fit)
        params = fit.param
        chisq = sum(fit.resid.^2)/2
    else
        function f(B, x)
            return B[1].*x .+ B[2]
        end
        linear = SciPy.odr.Model(f)
        fitdata = SciPy.odr.RealData(x_vals, y_vals, sx = x_err, sy = y_err)
        myodr = SciPy.odr.ODR(fitdata, linear, beta0=[1.0, 1.0])
        myoutput = myodr.run()
        cov = myoutput.cov_beta * myoutput.res_var
        params = myoutput.beta
        chisq = myoutput.res_var
    end
    @info "Fit Parameters", params
    @info "Reduced χ²", chisq
    range = LinRange(range[1], range[2], 1000)
    y = collect(range).*params[1] .+ params[2]
    if plotrange != :default
        range = plotrange
    end
    return (params=(m=params[1], c=params[2]), residuals=chisq, plot=plot!(range, y, ribbon = x -> sqrt(cov[4] + cov[1]*x^2 + 2*x*cov[2])))
end



function filter(data::DataFrame; β = :all, csw = :all)
    return DataFrames.filter([:β, :csw] => (β_v, csw_v) -> (csw_v == csw || csw == :all) && (β_v == β || β == :all), data)
end

function load_data(path::String)
    data = load_yaml_as_dataframe(path)
    data[!, :mpi2] = select(data, :g5_folded => ByRow(propagate_square) => :mpi2)[:, 1]
    data[!, :w0mpi] = select(data, [:w0, :g5_folded] => ByRow(propagate_product) => :w0mpi)[:, 1]
    data[!, :w0mpcac] = select(data, [:w0, :m_pcac] => ByRow(propagate_product) => :w0mpcac)[:, 1]
    data[!, :w0mpi2] = select(data, [:w0mpi] => ByRow(propagate_square) => :w0mpi2)[:, 1]
    data[!, :w0mv] = select(data, [:w0, :gk_folded] => ByRow(propagate_product) => :w0mv)[:, 1]
    data[!, :smearing_radius] = select(data, [:w0, :geometry] => ByRow((x, y) -> sqrt(8)*x/y[3]) => :smearing_radius)[:, 1]
    data[!, :L] = select(data, :geometry => ByRow(x -> last(x)) => :geometry)[:, 1]
    data[!, :mpiL] = select(data, [:L, :g5_folded] => ByRow((L, mpi) -> L*mpi[1]) => :mpiL)[:, 1]
    data[!, :T] = select(data, :geometry => ByRow(x -> first(x)) => :geometry)[:, 1]
    return data
end

