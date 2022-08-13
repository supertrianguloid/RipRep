using DrWatson
@quickactivate "RipRep"
using OffsetArrays
using DataFrames

@info "This is " * projectname() * " running from " * projectdir()

FLOAT_REGEX = r"([-+]?[0-9]*\.?[0-9]*)"
INTEGER_REGEX = r"([0-9]+)"

GLOBAL_SIMS = DataFrame(β=Float64[],
    Cˢʷ=Float64[],
    Mass=Float64[],
    L=Int64[],
    T=Int64[],
    SimulationType=String[],
    path=String[])

ensembledir() = datadir() * "/ensembles"
output_file_location(run_no) = ensembledir() * "/" * GLOBAL_SIMS[run_no, :path] * "/out_0"

function initialise()
    empty!(GLOBAL_SIMS)
    for sim in readdir(ensembledir())
        r = r"b" * FLOAT_REGEX * "_csw" * FLOAT_REGEX * "_m" * FLOAT_REGEX * "_L" * INTEGER_REGEX * "T" * INTEGER_REGEX * r"_(.+)"
        m = match(r, sim).captures
        push!(GLOBAL_SIMS, [map(parse, eltype.(eachcol(GLOBAL_SIMS))[1:end-1], m[1:end-1]); m[end]; sim])
    end
end

function list_ensembles()
    show(GLOBAL_SIMS[!, Not(:path)])
end

function load_outfile(run_no)
    return _convert_outformat_to_dataframe(read(output_file_location(run_no), String))
end

function load_processes(run_no)
    # Open the output file and split on new processes
    output = read(output_file_location(run_no), String)
    p = split(output, "[SYSTEM][0]Gauge group: SU(2)\n", keepempty=false)
    processes = Vector{DataFrame}(undef, length(p))
    # Populate a vector of DataFrames for each process
    for i in eachindex(p)
        processes[i] = _convert_outformat_to_dataframe(p[i])
    end
    return processes
end

function parse_output_file(run_no)

    processes = load_processes(run_no)

    global_metadata = extract_global_metadata(processes[1])

    process_metadata = Array{DataFrame}(undef, length(processes))

    for i in eachindex(processes)
        process_metadata[i] = process_metadata_to_dataframe(processes[i])
        process_metadata[i][1, :process_no] = i
    end

    process_metadata = vcat(process_metadata...)

    return global_metadata, process_metadata
end

function _check_integrator(p_meta::DataFrame)
    intp = dropmissing(select(p_meta, :process_no, names(p_meta, r"int_") .=> (x -> [i == 1 ? missing : x[i] != x[i-1] for i in axes(x, 1)]), renamecols=false))
    changes = intp[reduce(|, hcat([intp[:, i] .!= 0 for i in names(intp, r"int_")]...), dims=2)[:], :]
    return changes[:, :process_no]
end

function _convert_outformat_to_dataframe(process::AbstractString)
    run = DataFrame(name=String[], output=String[])
    line_regex = r"^\[([A-Z_]+)\]\[[0-9]+\](.*)"
    lines = split(process, "\n", keepempty=false)
    for line in lines
        m = match(line_regex, line).captures
        push!(run, m)
    end
    return run
end

function extract_global_metadata(first_process::DataFrame)
    ranlux = only(match(r"^RLXD (.+)", only(first_process[first_process.name.=="SETUP_RANDOM", :output])).captures)
    return DataFrame(ranlux=ranlux)
end


function process_metadata_to_dataframe(process::DataFrame)
    df = DataFrame(process_no=Int64[],
        warnings=Any[],
        int_l0_int=String[],
        int_l0_steps=Integer[],
        int_l1_int=String[],
        int_l1_steps=Integer[],
        int_l2_int=String[],
        int_l2_steps=Integer[]
    )
    # Helper function
    function extract_output(proc, name)
        p = proc[proc.name.==name, :output]
        return p == String[] ? nothing : p
    end

    # Extract warnings
    warnings = extract_output(process, "WARNING")

    # Extract integrator parameters
    integrator = []
    integrator_regex = r"Level " * INTEGER_REGEX * r": type = (o2mn|o4nm), steps = " * INTEGER_REGEX
    for i in extract_output(process, "INTEGRATOR")
        push!(integrator, match(integrator_regex, i).captures)
    end
    ints = Dict()
    steps = Dict()
    for i in integrator
        ints[i[1]] = String(i[2])
        steps[i[1]] = parse(Int, i[3])
    end

    push!(df, (process_no=0,
        warnings=warnings,
        int_l0_int=ints["0"],
        int_l0_steps=steps["0"],
        int_l1_int=ints["1"],
        int_l1_steps=steps["1"],
        int_l2_int=ints["2"],
        int_l2_steps=steps["2"],
    ))
    return df
end

initialise()


# function measure_acceptance(ensemble)
#     Nₐ = sum((ensemble.name .== "HMC") .& (ensemble.output .== "Configuration accepted."))
#     Nᵣ = sum((ensemble.name .== "HMC") .& (ensemble.output .== "Configuration rejected."))
#     return Nₐ/(Nₐ + Nᵣ)
# end

# function parse_process(process)

# end



# function parse_process(ensemble)
#     meas = DataFrame()
#     observables = Set()
#     for line in ensemble[ensemble.name .== "MAIN", :].output
#         m = match(r"^conf #[0-9]+ mass=" * fl * r" DEFAULT_SEMWALL TRIPLET ([a-zA-Z0-9_]+)=", line)
#         if typeof(m) != Nothing
#             push!(observables, m.captures[2])
#         end
#     end
#     return observables
# end


