using OffsetArrays
using DataFrames
using PrettyTables
using StatsBase

@info "RipRep parser activated."

FLOAT_REGEX_FULL = r"([+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?))"
FLOAT_REGEX = r"([-+]?[0-9]*\.?[0-9]*)"

INTEGER_REGEX = r"([0-9]+)"
TIME_REGEX = r"Trajectory #" * INTEGER_REGEX * r": generated in \[" * INTEGER_REGEX * " sec " * INTEGER_REGEX * r" usec\]"
CONF_REGEX = r"\[DeltaS = " * FLOAT_REGEX_FULL * r"\]\[" * r"exp\(-DS\) = " * FLOAT_REGEX_FULL * r"\]"
WF_REGEX = r"Configuration from .*n" * INTEGER_REGEX

runs_ignore_file(path::String) = path * "/RipRep/ignore_runs"
out_file(path::String) = path * "/out_0"
reports_dir() = "../reports"

mutable struct Ensemble
    global_metadata::DataFrame
    run_metadata::DataFrame
    data::DataFrame
    analysis::DataFrame
end

mutable struct WilsonFlow
    metadata::DataFrame
    data::DataFrame
    analysis::DataFrame
end

function load_runs(path::String)
    # Open the output file and split on new runs
    output = read(out_file(path), String)
    r = split(output, "[SYSTEM][0]Gauge group: SU(2)\n", keepempty=false)
    runs = Vector{DataFrame}(undef, length(r))
    # Populate a vector of DataFrames for each run
    for i in eachindex(r)
        runs[i] = _convert_outformat_to_dataframe(r[i])
    end
    return runs
end

function load_ensemble(path::String)
    g_meta, r_meta, data = parse_output_file(path)
    data[:, :gk] = (data[:, :g1] + data[:, :g2] + data[:, :g3])/3
    e = Ensemble(g_meta, r_meta, data, data)
    e.global_metadata.IntegratorChanges = [_check_integrator(e)]
    e.global_metadata.Acceptance = [_measure_acceptance(e)]
    e.global_metadata.DuplicatedConfigurations = [sort([i[1] for i in countmap(e.data[:, :conf_no]) if i[2] > 1])]
    name = split(path, "/")[end]
    r = r"b" * FLOAT_REGEX * "_csw" * FLOAT_REGEX * "_m" * FLOAT_REGEX * "_L" * INTEGER_REGEX * "T" * INTEGER_REGEX * r"_(.+)"
    m = match(r, name).captures
    e.global_metadata.β = [parse(Float64, m[1])]
    e.global_metadata.Cˢʷ = [parse(Float64, m[2])]
    e.global_metadata.Mass = [parse(Float64, m[3])]
    e.global_metadata.L = [parse(Int64, m[4])]
    e.global_metadata.T = [parse(Int64, m[5])]
    e.global_metadata.SimulationType = [string(m[6])]
    e.global_metadata = select(e.global_metadata, [:β, :Cˢʷ, :Mass, :L, :T, :SimulationType, :Confs, :Acceptance, :Runs, :IntegratorChanges, :MissingConfigurations, :DuplicatedConfigurations, :ranlux])
    @show e.global_metadata
    return e
end

function parse_output_file(path::String)

    runs = load_runs(path)
    runs = _prune_runs(runs)

    usable_runs = 1:length(runs)
    if(isfile(runs_ignore_file(path)))
        ignore = parse.(Int, split(strip(read(runs_ignore_file(path), String)), ','))
        usable_runs = []
        for i in 1:length(runs)
            if(i ∉ ignore)
                push!(usable_runs, i)
            end
        end
    end

    global_metadata = extract_global_metadata(first(runs[usable_runs]))

    run_metadata = DataFrame()

    for i in usable_runs
        push!(run_metadata, run_metadata_to_dataframe(runs[i])[1,:])
        run_metadata[end, :run_no] = i
    end

    data = Array{DataFrame}(undef, nrow(run_metadata))

    for i in 1:length(usable_runs)
        data[i] = extract_data_from_run(runs[usable_runs[i]])
        data[i][:, :run_no] .= usable_runs[i]
    end

    data = vcat(data...)

    global_metadata.Runs = [nrow(run_metadata)]
    global_metadata.Confs = [nrow(data)]
    global_metadata.MissingConfigurations = [data[([any(r) for r in eachrow(ismissing.(data))]), :].conf_no]

    return global_metadata, run_metadata, data
end

function _convert_outformat_to_dataframe(str::AbstractString)
    df = DataFrame(name=String[], output=String[])
    line_regex = r"^\[([a-zA-Z_]+)\]\[[0-9]+\](.*)"
    lines = split(str, "\n", keepempty=false)
    for line in lines
        try
            m = match(line_regex, line).captures
            push!(df, m)
        catch e
            @error "Error parsing line: " * line
        end
    end
    return df
end

function _prune_runs(runs::Vector{DataFrame})
    runs = runs[["ERROR" ∉ i.name for i in runs]]
    runs = runs[["Unknown integrator type" ∉ i.output for i in runs]]
    runs = runs[["Configuration accepted." ∈ i.output || "Configuration rejected." ∈ i.output for i in runs]]
    #runs = runs[["Process finalized." ∈ i.output for i in runs]]
    return runs
end

function extract_global_metadata(first_frame::DataFrame)
    ranlux = only(match(r"^RLXD (.+)", only(_extract_output(first_frame, "SETUP_RANDOM"))).captures)
    return DataFrame(ranlux=ranlux)
end

function _extract_output(run, name)
    p = run[run.name.==name, :output]
    return p == String[] ? nothing : p
end

function run_metadata_to_dataframe(run::DataFrame)
    df = DataFrame(run_no=Int64[],
        warnings=Any[],
        starting_conf=Union{Integer, Missing}[],
        no_confs=Union{Integer, Missing}[],
        integrator=String[]
    )

    # Extract warnings
    warnings = _extract_output(run, "WARNING")

    # Extract number of configurations
    starting_conf = missing
    ending_conf = missing
    try
        starting_conf = parse(Int, match(r"Trajectory #" * INTEGER_REGEX * "...", _extract_output(run, "MAIN") |> join).captures[1])
        ending_conf = parse(Int, match(r"Trajectory #" * INTEGER_REGEX * "...", _extract_output(run, "MAIN") |> reverse |> join).captures[1])
    catch e
        @error "Could not determine starting/ending confs for run"
        @show run
    end
    push!(df, (run_no=0,
        warnings=warnings,
        starting_conf = starting_conf,
        no_confs = ending_conf - starting_conf + 1,
        integrator=join(_extract_output(run, "INTEGRATOR"), "; ")
    ))

    # Extract configurations
    split(join(_extract_output(run, "MAIN"), '\n'), r"")
    return df
end

function extract_data_from_run(run::DataFrame)
    observables = [:g5, :g5_im, :id, :id_im, :g0, :g0_im, :g1, :g1_im, :g2, :g2_im, :g3, :g3_im, :g0g1, :g0g1_im, :g0g2, :g0g2_im, :g0g3, :g0g3_im, :g0g5, :g0g5_im, :g5g1, :g5g1_im, :g5g2, :g5g2_im, :g5g3, :g5g3_im, :g0g5g1, :g0g5g1_im, :g0g5g2, :g0g5g2_im, :g0g5g3, :g0g5g3_im, :g5_g0g5_re, :g5_g0g5_im]

    outer = []
    inner = []

    for i in run[:, :output]
        if(match(CONF_REGEX, i) !== nothing)
            push!(outer, inner)
            inner = []
        end
        push!(inner, i)
    end
    push!(outer, inner)

    confs = [join(i, '\n') for i in outer]

    df = DataFrame([Int64[], Int64[], Float64[], Union{Bool, Missing}[], Union{Float64, Missing}[], fill(Vector{Union{OffsetVector{Float64, Vector{Float64}}, Missing}}(), length(observables))...], [:conf_no, :run_no, :time, :accepted, :plaquette, observables...])
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

function _check_integrator(e::Ensemble)
    intp = dropmissing(select(e.run_metadata, :run_no, names(e.run_metadata, r"integrator") .=> (x -> [i == 1 ? missing : x[i] != x[i-1] for i in axes(x, 1)]), renamecols=false))
    changes = intp[reduce(|, hcat([intp[:, i] .!= 0 for i in names(intp, r"integrator")]...), dims=2)[:], :]
    return _run_to_first_conf(e, changes[:, :run_no])
end

function _run_to_first_conf(e::Ensemble, run_no::Vector{Int64})
    confs = []
    for i in run_no
        try
            push!(confs, e.data[e.data.run_no .== i, :][1, :conf_no])
        catch e
            @error "Error finding first configuration of run " * string(i)
        end
    end
    return confs
end

function _run_to_first_conf(e::Ensemble, run_no::Int64)
    return _run_to_first_conf(e, [run_no])
end

function _measure_acceptance(e::Ensemble)
    data = e.data
    if(only(e.global_metadata.IntegratorChanges) != [])
        data = e.data[maximum(only(e.global_metadata.IntegratorChanges)):end, :]
    end
    return sum(data[:, :accepted])/nrow(data)
end

function load_wilson_flow(path)
    output = split(read(path, String), '\n')
    
    data = DataFrame([Int64[], Float64[], Float64[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]], [:conf_no, :time, :plaquette, :t, :E, :t²E, :Esym, :t²Esym, :TC])

    outer = []
    inner = []

    for i in output
        if(match(WF_REGEX, i) !== nothing)
            push!(outer, inner)
            inner = []
        end
        push!(inner, i)
    end
    push!(outer, inner)
    outer = join.(outer, '\n')

    outer = _convert_outformat_to_dataframe.(outer)

    metadata = _extract_wf_metadata(outer[1])

    for conf in outer[2:end]
        conf_no = parse(Int, only(match(WF_REGEX, only(_extract_output(conf, "MAIN"))).captures))
        time = match(r"done \[" * INTEGER_REGEX * " sec " * INTEGER_REGEX * r" usec\]", only(_extract_output(outer[2], "TIMING"))).captures
        time = parse(Int, time[1]) + 1e-6*parse(Int, time[2])
        plaquette = parse(Float64, match("Plaquette="*FLOAT_REGEX_FULL, _extract_output(conf, "IO")[2]).captures[1])
        wf = []
        for flow in _extract_output(conf, "WILSONFLOW")
            a = match(r"WF \(t,E,t2\*E,Esym,t2\*Esym,TC\) = (.*)", flow)
            push!(wf, parse.(Float64, split(a.captures[1], ' ')))
        end
        wf = hcat(wf...)'
        push!(data, (conf_no = conf_no,
                    time = time,
                    plaquette = plaquette,
                    t = wf[:, 1],
                    E = wf[:, 2],
                    t²E = wf[:, 3],
                    Esym = wf[:, 4],
                    t²Esym = wf[:, 5], 
                    TC = wf[:, 6]
                    ))
    end                 
    return WilsonFlow(metadata, data, data)
end

function _extract_wf_metadata(first_frame::DataFrame)
    ranlux = _only_match(r"^RLXD (.+)$", _extract_output(first_frame, "SETUP_RANDOM"))
    integrator = _only_match(r"WF integrator: (.+)$", _extract_output(first_frame, "MAIN"))
    tmax = parse(Float64, _only_match(r"tmax: (.+)$", _extract_output(first_frame, "MAIN")))
    N_meas = parse(Int64, _only_match(r"WF number of measures: (.+)$", _extract_output(first_frame, "MAIN")))
    ϵ₀ = parse(Float64, _only_match(r"WF initial epsilon: (.+)$", _extract_output(first_frame, "MAIN")))
    Δ = parse(Float64, _only_match(r"WF delta: (.+)$", _extract_output(first_frame, "MAIN")))
    dt = parse(Float64, _only_match(r"WF measurement interval dt : (.+)$", _extract_output(first_frame, "MAIN")))

    return DataFrame(integrator=integrator,
                    ranlux=ranlux,
                    tmax=tmax,
                    N_meas=N_meas,
                    ϵ₀=ϵ₀,
                    Δ=Δ,
                    dt=dt)
end

