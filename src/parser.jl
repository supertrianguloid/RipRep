using OffsetArrays
using DataFrames
using PrettyTables
using StatsBase

include("utilities.jl")

FLOAT_REGEX_FULL = r"([+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?))"
FLOAT_REGEX = r"([-+]?[0-9]*\.?[0-9]*)"

INTEGER_REGEX = r"([0-9]+)"
TIME_REGEX = r"Trajectory #" * INTEGER_REGEX * r": generated in \[" * INTEGER_REGEX * " sec " * INTEGER_REGEX * r" usec\]"
CONF_REGEX = r"\[DeltaS = " * FLOAT_REGEX_FULL * r"\]\[" * r"exp\(-DS\) = " * FLOAT_REGEX_FULL * r"\]"
WF_REGEX = r"Configuration from .*n" * INTEGER_REGEX

runs_ignore_file(path::String) = path * "/RipRep/ignore_runs"
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
    output = read(path, String)
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
    data[:, :dg5_g0g5_re] = _corr_derivative(data[:, :g5_g0g5_re])
    data[:, :gk_folded] = _fold(data[:, :gk])
    data[:, :g5_folded] = _fold(data[:, :g5])
    data[:, :id_folded] = _fold(data[:, :id])
    data[:, :g5_g0g5_re_folded] = _fold(data[:, :g5_g0g5_re], symm = false)
    data[:, :dg5_g0g5_re_folded] = _corr_derivative(data[:, :g5_g0g5_re_folded])

    e = Ensemble(g_meta, r_meta, data, data)
    e.global_metadata.IntegratorChanges = [_check_integrator(e)]
    e.global_metadata.Acceptance = [_measure_acceptance(e)]
    e.global_metadata.DuplicatedConfigurations = [sort([i[1] for i in countmap(e.data[:, :conf_no]) if i[2] > 1])]
    
    e.global_metadata.SimulationType = ["expclv"]
    e.global_metadata = select(e.global_metadata, [:β, :Cˢʷ, :Mass, :L, :T, :SimulationType, :Confs, :Acceptance, :Runs, :IntegratorChanges, :MissingConfigurations, :DuplicatedConfigurations, :ranlux])
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
    if(all(data.eig_max .=== missing) && all(data.eig_min .=== missing))
        select!(data, Not([:eig_max, :eig_min]))
    end
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
    ranlux = missing
    try
        ranlux = only(match(r"^RLXD (.+)", only(_extract_output(first_frame, "SETUP_RANDOM"))).captures)
    catch e
        print("Error extracing SETUP_RANDOM")
    end
    action = _extract_output(first_frame, "ACTION")
    T = parse(Int64, split(split(_extract_output(first_frame,"GEOMETRY_INIT")[1])[end], "x")[1])
    L = parse(Int64, split(split(_extract_output(first_frame,"GEOMETRY_INIT")[1])[end], "x")[end])
    β = parse(Float64, match(r"beta = ([\+\-0-9.]+)", action[2]).captures[1])
    Mass = parse(Float64, match(r"mass = ([\+\-0-9.]+)", action[3]).captures[1])
    Cˢʷ = parse(Float64, only(match(r"reset to csw = ([\+\-0-9.]+)", _extract_output(first_frame, "CLOVER")[end]).captures))
    return DataFrame(ranlux=ranlux, L=L, T=T, β=β, Mass=Mass, Cˢʷ=Cˢʷ)
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

    df = DataFrame([Int64[], Int64[], Float64[], Union{Bool, Missing}[], Union{Float64, Missing}[], Union{Float64, Missing}[],Union{Float64, Missing}[],Union{Float64, Missing}[],Union{Float64, Missing}[], Union{Float64, Missing}[],Union{Float64, Missing}[],Union{Float64, Missing}[], fill(Vector{Union{OffsetVector{Float64, Vector{Float64}}, Missing}}(), length(observables))...], [:conf_no, :run_no, :time, :accepted, :dS, :Polyakov_0, :Polyakov_1, :Polyakov_2, :Polyakov_3, :eig_min, :eig_max, :plaquette, observables...])
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
            push!(d, parse(Float64, only(match(r"Plaquette: " * r"(.*)", i).captures)))
        catch e
            push!(d, missing)
        end

        status = missing
        try
            status = (only(match(r"Configuration (rejected|accepted).", i).captures) == "accepted")
        catch e
        end

        dS = missing
        try
            dS = only(match(r"\[DeltaS = (.*)\]\[exp", i).captures)
            dS = parse(Float64, dS)
        catch e
        end

        Polyakov_0 = missing
        try
            Polyakov_0 = only(match(r"Polyakov direction 0 = (.*)", i).captures)
            Polyakov_0 = parse(Float64, split(Polyakov_0)[1])
        catch e
        end
        Polyakov_1 = missing
        try
            Polyakov_1 = only(match(r"Polyakov direction 1 = (.*)", i).captures)
            Polyakov_1 = parse(Float64, split(Polyakov_1)[1])
        catch e
        end
        Polyakov_2 = missing
        try
            Polyakov_2 = only(match(r"Polyakov direction 2 = (.*)", i).captures)
            Polyakov_2 = parse(Float64, split(Polyakov_2)[1])
        catch e
        end
        Polyakov_3 = missing
        try
            Polyakov_3 = only(match(r"Polyakov direction 3 = (.*)", i).captures)
            Polyakov_3 = parse(Float64, split(Polyakov_3)[1])
        catch e
        end

        eig_min = missing
        try
            eig_min = only(match(r"Eig 0 = (.*)", i).captures)
            eig_min = parse(Float64, eig_min)
        catch e
        end

        eig_max = missing
        try
            eig_max = only(match(r"Max_eig = (.*)", i).captures)
            eig_max = parse(Float64, split(eig_max)[1])
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
        push!(df, [confno, 0, time, status, dS, Polyakov_0, Polyakov_1, Polyakov_2, Polyakov_3, eig_min, eig_max, d...])
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

    outer = outer[["TIMING" ∈ i.name for i in outer]]


    for conf in outer
        conf_no = parse(Int, only(match(WF_REGEX, only(_extract_output(conf, "MAIN"))).captures))
        time = match(r"done \[" * INTEGER_REGEX * " sec " * INTEGER_REGEX * r" usec\]", only(_extract_output(conf, "TIMING"))).captures
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
    wf = WilsonFlow(metadata, data, data)
    W = Vector{Vector{Float64}}()
    for i in wf.data.t²E
        w = Vector{Float64}()
        for j in 2:(length(i)-1)
            push!(w, (i[j+1] - i[j-1]))
        end
        w = w .* (wf.data[1, :t][2:(length(i)-1)])
        push!(W, w)
    end
    W = W./(2*wf.metadata.dt)
    wf.analysis[:, :W] = wf.data[:, :W] = W
    return wf
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
function load_sf(path)
    output = split(read(path, String), '\n')

    outer = []
    inner = []

    for i in output
        if(match(CONF_REGEX, i) !== nothing)
            push!(outer, inner)
            inner = []
        end
        push!(inner, i)
    end
    push!(outer, inner)
    outer = join.(outer, '\n')

    outer = _convert_outformat_to_dataframe.(outer)
    df = DataFrame(confno = Int[], accepted = Bool[], time = Union{Float64, Missing}[], dH =  Union{Float64, Missing}[], e⁻ᵈᴴ = Union{Float64, Missing}[], plaquette = Union{Float64, Missing}[])

    for i in outer[2:end]
        dH = parse(Float64, split(split(_extract_output(i, "HMC")[1])[3], "]")[1])
        e⁻ᵈᴴ = parse(Float64, split(split(_extract_output(i, "HMC")[1])[5], "]")[1])
        accepted = (only(match(r"Configuration (rejected|accepted).", _extract_output(i, "HMC")[2]).captures) == "accepted")
        plaq = missing
        try
            plaq = parse(Float64, split(only([i for i in _extract_output(i, "MAIN") if occursin("Plaquette:", i)]))[2])
        catch e
        end
        t = match(TIME_REGEX, _extract_output(i, "MAIN")[1]).captures
        time = parse(Int, t[2]) + 1e-6*parse(Int, t[3])
        confno = parse(Int, t[1])
#        push!(df, (confno, time, parse(Float64, hmc[1]), parse(Float64, hmc[2])))
        push!(df, (confno, accepted, time, dH, e⁻ᵈᴴ, plaq))
    end
    return df
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

function _corr_derivative(correlator; h = 1)
    corr = []
    for c in correlator
        if c === missing
            push!(corr, missing)
            continue
        end
        push!(corr, d(c))
    end
    return corr
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
