using DrWatson
@quickactivate "RipRep"
using OffsetArrays
using DataFrames
using PrettyTables

@info "This is " * projectname() * " running from " * projectdir()

FLOAT_REGEX_FULL = r"([+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?))"
FLOAT_REGEX = r"([-+]?[0-9]*\.?[0-9]*)"

INTEGER_REGEX = r"([0-9]+)"
TIME_REGEX = r"Trajectory #" * INTEGER_REGEX * r": generated in \[" * INTEGER_REGEX * " sec " * INTEGER_REGEX * r" usec\]"
CONF_REGEX = r"\[DeltaS = " * FLOAT_REGEX_FULL * r"\]\[" * r"exp\(-DS\) = " * FLOAT_REGEX_FULL * r"\]"
WF_REGEX = r"Configuration from .*n" * INTEGER_REGEX

GLOBAL_SIMS = DataFrame(β=Float64[],
    Cˢʷ=Float64[],
    Mass=Float64[],
    L=Int64[],
    T=Int64[],
    SimulationType=String[],
    path=String[])

GLOBAL_WF = DataFrame(β=Float64[],
    Cˢʷ=Float64[],
    Mass=Float64[],
    L=Int64[],
    T=Int64[],
    SimulationType=String[],
    path=String[])

struct Ensemble
    global_metadata::DataFrame
    run_metadata::DataFrame
    data::DataFrame
end

struct WilsonFlow
    metadata::DataFrame
    data::DataFrame
end

ensembledir() = datadir() * "/ensembles"
reportsdir() = projectdir() * "/reports"

output_file_location(ens_no) = ensembledir() * "/" * GLOBAL_SIMS[ens_no, :path] * "/out_0"

wfdir() = datadir() * "/WF"
wf_file_location(run_no) = wfdir() * "/" * GLOBAL_WF[run_no, :path] * "/wilsonflow.out"

function initialise()
    empty!(GLOBAL_SIMS)
    for sim in readdir(ensembledir())
        r = r"b" * FLOAT_REGEX * "_csw" * FLOAT_REGEX * "_m" * FLOAT_REGEX * "_L" * INTEGER_REGEX * "T" * INTEGER_REGEX * r"_(.+)"
        m = match(r, sim).captures
        push!(GLOBAL_SIMS, [map(parse, eltype.(eachcol(GLOBAL_SIMS))[1:end-1], m[1:end-1]); m[end]; sim])
    end
    empty!(GLOBAL_WF)
    for sim in readdir(wfdir())
        r = r"WF_b" * FLOAT_REGEX * "_csw" * FLOAT_REGEX * "_m" * FLOAT_REGEX * "_L" * INTEGER_REGEX * "T" * INTEGER_REGEX * r"_(.+)"
        m = match(r, sim).captures
        push!(GLOBAL_WF, [map(parse, eltype.(eachcol(GLOBAL_WF))[1:end-1], m[1:end-1]); m[end]; sim])
    end
end

function list_ensembles()
    return GLOBAL_SIMS[!, Not(:path)]
end

function list_wilsonflows()
    return GLOBAL_WF[!, Not(:path)]
end

function load_ensemble(ens_no)
    return _convert_outformat_to_dataframe(read(output_file_location(ens_no), String))
end

function load_runs(ens_no)
    # Open the output file and split on new runs
    output = read(output_file_location(ens_no), String)
    r = split(output, "[SYSTEM][0]Gauge group: SU(2)\n", keepempty=false)
    runs = Vector{DataFrame}(undef, length(r))
    # Populate a vector of DataFrames for each run
    for i in eachindex(r)
        runs[i] = _convert_outformat_to_dataframe(r[i])
    end
    return runs
end

function load_wilson_flow(wf_no)
    output = split(read(wf_file_location(wf_no), String), '\n')
    
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
    return WilsonFlow(metadata, data)
end

function parse_output_file(ensemble_no)

    runs = load_runs(ensemble_no)
    runs = runs[["ERROR" ∉ i.name for i in runs]]

    global_metadata = extract_global_metadata(runs[1])

    run_metadata = Array{DataFrame}(undef, length(runs))

    for i in eachindex(runs)
        run_metadata[i] = run_metadata_to_dataframe(runs[i])
        run_metadata[i][1, :run_no] = i
    end

    run_metadata = vcat(run_metadata...)

    data = Array{DataFrame}(undef, nrow(run_metadata))

    for i in eachindex(runs)
        data[i] = extract_data_from_run(runs[i])
        data[i][:, :run_no] .= i
    end

    data = vcat(data...)

    global_metadata.Runs = [nrow(run_metadata)]
    global_metadata.Confs = [nrow(data)]
    global_metadata.IntegratorChanges = [_check_integrator(run_metadata)]
    global_metadata.MissingConfigurations = [data[([any(r) for r in eachrow(ismissing.(data))]), :].conf_no]
    global_metadata.Acceptance = [_measure_acceptance(data)]

    return global_metadata, run_metadata, data
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
    intp = dropmissing(select(p_meta, :run_no, names(p_meta, r"int_") .=> (x -> [i == 1 ? missing : x[i] != x[i-1] for i in axes(x, 1)]), renamecols=false))
    changes = intp[reduce(|, hcat([intp[:, i] .!= 0 for i in names(intp, r"int_")]...), dims=2)[:], :]
    return changes[:, :run_no]
end

function _convert_outformat_to_dataframe(str::AbstractString)
    df = DataFrame(name=String[], output=String[])
    line_regex = r"^\[([A-Z_]+)\]\[[0-9]+\](.*)"
    lines = split(str, "\n", keepempty=false)
    for line in lines
        m = match(line_regex, line).captures
        push!(df, m)
    end
    return df
end

function _extract_output(proc, name)
    p = proc[proc.name.==name, :output]
    return p == String[] ? nothing : p
end

function extract_global_metadata(first_frame::DataFrame)
    ranlux = only(match(r"^RLXD (.+)", only(_extract_output(first_frame, "SETUP_RANDOM"))).captures)
    return DataFrame(ranlux=ranlux)
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


function run_metadata_to_dataframe(run::DataFrame)
    df = DataFrame(run_no=Int64[],
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
    warnings = _extract_output(run, "WARNING")

    # Extract integrator parameters
    integrator = []
    integrator_regex = r"Level " * INTEGER_REGEX * r": type = (o2mn|o4mn), steps = " * INTEGER_REGEX
    for i in _extract_output(run, "INTEGRATOR")
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
        starting_conf = parse(Int, match(r"Trajectory #" * INTEGER_REGEX * "...", _extract_output(run, "MAIN") |> join).captures[1])
        ending_conf = parse(Int, match(r"Trajectory #" * INTEGER_REGEX * "...", _extract_output(run, "MAIN") |> reverse |> join).captures[1])
    catch e
        print(run)
    end
    push!(df, (run_no=0,
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
    split(join(_extract_output(run, "MAIN"), '\n'), r"")
    return df
end

function load_ensemble(ensemble_no)
    g_meta, r_meta, data = parse_output_file(ensemble_no)

    g_meta = hcat(DataFrame(GLOBAL_SIMS[ensemble_no, :]), g_meta)

    return Ensemble(g_meta, r_meta, data)
end

function _ensemble_to_latex_string(ens)
    return "\$\\beta = " * string(only(ens.global_metadata.β)) * ",\\ C_{SW} = " * string(only(ens.global_metadata.Cˢʷ)) * ",\\ m = " * string(only(ens.global_metadata.Mass)) * ",\\ V = " * string(only(ens.global_metadata.T)) * "\\times" * string(only(ens.global_metadata.L)) * "^3\$"
end

function _measure_acceptance(data)
    return sum(data[:, :accepted])/nrow(data)
end

function summarise_ensembles(method = :show)
    df = DataFrame()
    for i in 1:nrow(GLOBAL_SIMS)
        push!(df, load_ensemble(i).global_metadata[1, :])
    end
    if(method == :save)
        open(reportsdir() * "/ensembles.md", "w") do f
            pretty_table(f,df, tf=PrettyTables.tf_markdown, nosubheader=true)
        end
    end
    return df
end

function _only_match(regex, output)
    m = match.(regex, output)
    return only(only(m[m .!= nothing]).captures)
end

initialise()