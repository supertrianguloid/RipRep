using DrWatson
@quickactivate "RipRep"
using OffsetArrays
using DataFrames

@info "This is " * projectname() * " running from " * projectdir()

function parse_sims()
    sims = DataFrame(β = Float64[], Cˢʷ = Float64[], Mass = Float64[], L = Int64[], T = Int64[], SimulationType = String[])
    for sim in readdir(datadir()*"/ensembles")
        println(sim)
        fl = r"([-+]?[0-9]*\.?[0-9]*)"
        int = r"([0-9]+)"
        r = r"b"*fl*"_csw"*fl*"_m"*fl*"_L"*int*"T"*int*r"_(.+)"
        m = match(r, sim).captures

        push!(sims, [map(parse, eltype.(eachcol(sims))[1:end-1], m[1:end-1]); m[end]])
    end
    print(sims)
end

parse_sims()
