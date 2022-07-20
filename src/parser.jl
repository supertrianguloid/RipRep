using DrWatson
@quickactivate "RipRep"
using OffsetArrays
using DataFrames

@info "This is " * projectname() * " running from " * projectdir()

function initialise()
    sims = DataFrame(β = Float64[], Cˢʷ = Float64[], Mass = Float64[], L = Int64[], T = Int64[], SimulationType = String[], path = String[])
    for sim in readdir(datadir()*"/ensembles")
        fl = r"([-+]?[0-9]*\.?[0-9]*)"
        int = r"([0-9]+)"
        r = r"b"*fl*"_csw"*fl*"_m"*fl*"_L"*int*"T"*int*r"_(.+)"
        m = match(r, sim).captures
        push!(sims, [map(parse, eltype.(eachcol(sims))[1:end-1], m[1:end-1]); m[end]; sim])
    end
    return sims
end


function load_simulation(sims, runno)
    run = DataFrame(name = String[], output = String[])
    line_regex = r"^\[([A-Z_]+)\]\[[0-9]+\](.*)"
    f = open(datadir()*"/ensembles/"*sims[runno, :path]*"/out_0")
    for line in eachline(f)
        m = match(line_regex, line).captures
        push!(run, m)
    end
    return run    
end

function measure_acceptance(ensemble)
    Nₐ = sum((ensemble.name .== "HMC") .& (ensemble.output .== "Configuration accepted."))
    Nᵣ = sum((ensemble.name .== "HMC") .& (ensemble.output .== "Configuration rejected."))
    return Nₐ/(Nₐ + Nᵣ)
end
