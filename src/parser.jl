using DrWatson
@quickactivate "RipRep"
using OffsetArrays
using DataFrames

@info "This is " * projectname() * " running from " * projectdir()

FLOAT_REGEX_FULL = r"([+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?))"
FLOAT_REGEX = r"([-+]?[0-9]*\.?[0-9]*)"

INTEGER_REGEX = r"([0-9]+)"
TIME_REGEX = r"Trajectory #" * INTEGER_REGEX * r": generated in \[" * INTEGER_REGEX * " sec " * INTEGER_REGEX * r" usec\]"
CONF_REGEX = r"\[DeltaS = " * FLOAT_REGEX_FULL * r"\]\[" * r"exp\(-DS\) = " * FLOAT_REGEX_FULL * r"\]"

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

    data = Array{DataFrame}(undef, nrow(process_metadata))

    for i in eachindex(processes)
        data[i] = extract_data_from_process(processes[i])
        data[i][:, :proc_no] .= i
    end

    data = vcat(data...)

    # for i in 1:nrow(data)
    #     data[i, :conf_no] = i
    # end

    md = data[([any(r) for r in eachrow(ismissing.(data))]), :].conf_no


    if md != []
        @info "Missing data in configuration $md"
    end

    return global_metadata, process_metadata, data
end

function extract_data_from_process(process::DataFrame)
    observables = [:g5, :g5_im, :id, :id_im, :g0, :g0_im, :g1, :g1_im, :g2, :g2_im, :g3, :g3_im, :g0g1, :g0g1_im, :g0g2, :g0g2_im, :g0g3, :g0g3_im, :g0g5, :g0g5_im, :g5g1, :g5g1_im, :g5g2, :g5g2_im, :g5g3, :g5g3_im, :g0g5g1, :g0g5g1_im, :g0g5g2, :g0g5g2_im, :g0g5g3, :g0g5g3_im, :g5_g0g5_re, :g5_g0g5_im]

    outer = []
    inner = []

    for i in process[:, :output]
        if(match(CONF_REGEX, i) !== nothing)
            push!(outer, inner)
            inner = []
        end
        push!(inner, i)
    end
    push!(outer, inner)

    confs = [join(i, '\n') for i in outer]

    df = DataFrame([Int64[], Int64[], Float64[], Union{Bool, Missing}[], Union{Float64, Missing}[], fill(Vector{Union{OffsetVector{Float64, Vector{Float64}}, Missing}}(), length(observables))...], [:conf_no, :proc_no, :time, :accepted, :plaquette, observables...])
    for i in confs[2:end]

        time = missing
        confno = missing
        try
            t = match(TIME_REGEX, i).captures
            time = parse(Int, t[2]) + 1e-6*parse(Int, t[3])
            confno = parse(Int, t[1])
        catch e
            println(i)
        end

        d = []
        try
            push!(d, parse(Float64, only(match(r"Plaquette: " * FLOAT_REGEX, i).captures)))
        catch e
            push!(d, missing)
        end

        status = missing
        try
            status = (only(match(r"Configuration (rejected|accepted).", i).captures) == "accepted")
        catch e
        end

        for obs in observables
            try
                # v = parse.(Float64, split(only(match(r"conf #" * String(confno) * r".*" * String(obs) * r"=(.*)", i).captures), ' ', keepempty=false))
                v = parse.(Float64, split(only(match(String(obs) * r"=(.*)", i).captures), ' ', keepempty=false))
                push!(d, OffsetArray(v, 0:length(v) - 1))
                
            catch e
                push!(d, missing)
            end
        end
        push!(df, [confno, 0, time, status, d...])
    end
    return df
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

function _extract_output(proc, name)
    p = proc[proc.name.==name, :output]
    return p == String[] ? nothing : p
end

function extract_global_metadata(first_process::DataFrame)
    ranlux = only(match(r"^RLXD (.+)", only(first_process[first_process.name.=="SETUP_RANDOM", :output])).captures)
    return DataFrame(ranlux=ranlux)
end


function process_metadata_to_dataframe(process::DataFrame)
    df = DataFrame(process_no=Int64[],
        warnings=Any[],
        starting_conf=Union{Integer, Missing}[],
        no_confs=Union{Integer, Missing}[],
        int_l0_int=String[],
        int_l0_steps=Integer[],
        int_l1_int=String[],
        int_l1_steps=Integer[],
        int_l2_int=String[],
        int_l2_steps=Integer[]
    )

    # Extract warnings
    warnings = _extract_output(process, "WARNING")

    # Extract integrator parameters
    integrator = []
    integrator_regex = r"Level " * INTEGER_REGEX * r": type = (o2mn|o4mn), steps = " * INTEGER_REGEX
    for i in _extract_output(process, "INTEGRATOR")
        try
            push!(integrator, match(integrator_regex, i).captures)
        catch e
            @info i
        end
    end
    ints = Dict()
    steps = Dict()
    for i in integrator
        ints[i[1]] = String(i[2])
        steps[i[1]] = parse(Int, i[3])
    end

    # Extract number of configurations
    starting_conf = missing
    ending_conf = missing
    try
        starting_conf = parse(Int, match(r"Trajectory #" * INTEGER_REGEX * "...", _extract_output(process, "MAIN") |> join).captures[1])
        ending_conf = parse(Int, match(r"Trajectory #" * INTEGER_REGEX * "...", _extract_output(process, "MAIN") |> reverse |> join).captures[1])
    catch e
        print(process)
    end
    push!(df, (process_no=0,
        warnings=warnings,
        starting_conf = starting_conf,
        no_confs = ending_conf - starting_conf + 1,
        int_l0_int=ints["0"],
        int_l0_steps=steps["0"],
        int_l1_int=ints["1"],
        int_l1_steps=steps["1"],
        int_l2_int=ints["2"],
        int_l2_steps=steps["2"],
    ))

    # Extract configurations
    split(join(_extract_output(process, "MAIN"), '\n'), r"")
    return df
end

function load_run(run_no)
    g_meta, p_meta, data = parse_output_file(run_no)
    print("\n")
    println("Number of processes: ", nrow(p_meta))
    println("Total number of configurations: ", nrow(data))
    println("Integrator parameters change in process numbers: ", _check_integrator(p_meta))
    return data
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


