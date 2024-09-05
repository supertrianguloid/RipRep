using DataFrames
using OffsetArrays
using SHA

include("utilities.jl")

RUN_BEGIN = "Gauge group: SU(2)"
TRAJ_BEGIN_REGEX = r"Trajectory #(.*)\.\.\.$"
TRAJ_GENERATED_REGEX = r"Trajectory #(.*): generated in \[(.*) sec (.*) usec\]$"
ACCEPT_REJECT_REGEX = r"^Configuration (rejected|accepted).$"
DS_REGEX = r"^\[DeltaS = (.*)\]\[exp\(-DS\) = (.*)\]$"
MEAS_REGEX = r"^conf #(.*) mass=(.*) DEFAULT_SEMWALL TRIPLET (.*)= (.*) $"
WF_TRAJ_BEGIN_REGEX = r"^Configuration from (.*)$"
WF_TRAJ_BEGIN_REGEX_RUN_NUMBER = r"^Configuration from .*n([0-9]+)$"
WF_MEASUREMENT_REGEX = r"^WF \(t,E,t2\*E,Esym,t2\*Esym,TC\) = (.*)$"
WF_TIMING_REGEX = r".* done \[(.*) sec (.*) usec\]$"
MVM_REGEX = r"Trajectory .* ; (.*)$"

CORRELATORS = [:g5, :g5_im, :id, :id_im, :g0, :g0_im, :g1, :g1_im, :g2, :g2_im, :g3, :g3_im, :g0g1, :g0g1_im, :g0g2, :g0g2_im, :g0g3, :g0g3_im, :g0g5, :g0g5_im, :g5g1, :g5g1_im, :g5g2, :g5g2_im, :g5g3, :g5g3_im, :g0g5g1, :g0g5g1_im, :g0g5g2, :g0g5g2_im, :g0g5g3, :g0g5g3_im, :g5_g0g5_re, :g5_g0g5_im]

mutable struct Ensemble
    global_metadata::Dict
    run_metadata::DataFrame
    data::DataFrame
    analysis::DataFrame
end

mutable struct WilsonFlow
    metadata::Dict
    data::DataFrame
    analysis::DataFrame
end

function load_output_file_as_dataframe(path::String)
    df = DataFrame(lineno=Int[], name=String[], loglevel=Int[], output=String[])
    line_regex = r"^\[([a-zA-Z_ ]+)\]\[(-?[0-9]+)\](.*)"
    lineno = 0
    emptylines = 0
    for line in eachline(open(path))
        if line != ""
	    lineno += 1
            try
                m = match(line_regex, line).captures
                push!(df, [lineno, m[1], parse(Int, m[2]), m[3]])
            catch e
		@error "Error parsing output file on line: " * (line+emptylines)
                @show lineno
            end
	else
            emptylines += 1
        end
    end
    return df
end

function extract_global_metadata(path::String, output_df::DataFrame, data::DataFrame, run_metadata::DataFrame)
    global_metadata = Dict{Any, Any}(:rng => missing)

    global_metadata[:nconfs] = length(data.confno)

    global_metadata[:csw] = only(union(run_metadata.csw))
    global_metadata[:β] = only(union(run_metadata[:, :β]))
    global_metadata[:m0] = only(union(run_metadata[:, :m0]))
    try
        global_metadata[:mu] = last(union(run_metadata[:, :mu]))
    catch e
    end

    runs_where_the_integrator_changes = filter(:integrator => integrator -> integrator == true, dropmissing(select(run_metadata, :run_number, names(run_metadata, :integrator) .=>  (x -> [i == 1 ? missing : x[i] != x[i-1] for i in axes(x, 1)]), renamecols=false))).run_number

    runs_where_the_rng_is_reseeded = filter(:rng => rng -> rng != [], dropmissing(select(run_metadata, :run_number, names(run_metadata, :rng)))).run_number

    global_metadata[:reseeded_confs] = [run_to_first_conf(data, i) for i in runs_where_the_rng_is_reseeded]

    if length(global_metadata[:reseeded_confs]) > 1
        @warn "RNG reseeded until conf " * string(last(global_metadata[:reseeded_confs]))
    end
    
    global_metadata[:integrator_changes] = [run_to_first_conf(data, i) for i in runs_where_the_integrator_changes]
                                                
    global_metadata[:geometry] = only(union(run_metadata.geometry))

    global_metadata[:path] = path

    nan_confs = find_nans(data)
    
    if !isempty(nan_confs)
        @warn "The run contains NaNs"
    end

    global_metadata[:nan_confs] = nan_confs

    open(path) do f
        global_metadata[:sha256] = bytes2hex(sha2_256(f))
    end
    
    return global_metadata
end

function split_output_dataframe_into_runs(output_df::DataFrame; keep_runs_without_trajectories=false)
    runs = []
    nruns = output_df[end, :run_number]
    for run_number in 1:nruns
        run = filter(:run_number => x -> x == run_number, output_df, view=true)
        ntraj = nrow(filter(:output => x -> !isnothing(match(TRAJ_BEGIN_REGEX, x)), run, view=true))
        if keep_runs_without_trajectories || ntraj > 0
            push!(runs, run)
        end
    end
    return runs 
end

# TODO Needs to be able to handle runs with < 3 trajectories
function split_run_dataframe_into_trajectories(runs)
    trajectories = []
    for run_number in eachindex(runs)
        run_df = runs[run_number]
        traj_boundaries = filter([:name, :output] => (name, output) -> name == "MAIN" && !isnothing(match(TRAJ_BEGIN_REGEX, output)), run_df, view=true).lineno
        if length(traj_boundaries) == 0
            continue
        end
        if length(traj_boundaries) == 1
            push!(trajectories, run_df)
            continue
        end
        first_traj = only(describe(run_df, cols=:lineno, :min).min):(traj_boundaries[2] - 1)
        trajectory = filter(:lineno => lineno -> lineno ∈ first_traj, run_df, view=true)
        push!(trajectories, trajectory)
        for i in 2:(length(traj_boundaries) - 1)
            trajectory = filter(:lineno => lineno -> lineno ∈ traj_boundaries[i]:(traj_boundaries[i+1] - 1), run_df, view=true)
            push!(trajectories, trajectory)
        end
        last_traj = last(traj_boundaries):only(describe(run_df, cols=:lineno, :max).max)
        trajectory = filter(:lineno => lineno -> lineno ∈ last_traj, run_df, view=true)
        push!(trajectories, trajectory)
    end
    return trajectories
end

function split_wf_dataframe_into_trajectories(runs)
    trajectories = []
    for run_number in eachindex(runs)
        run_df = runs[run_number]
        traj_boundaries = filter([:name, :output] => (name, output) -> name == "MAIN" && !isnothing(match(WF_TRAJ_BEGIN_REGEX, output)), run_df, view=true).lineno
        if length(traj_boundaries) == 1
            push!(trajectories, run_df)
            continue
        end
        first_traj = only(describe(run_df, cols=:lineno, :min).min):(traj_boundaries[2] - 1)
        trajectory = filter(:lineno => lineno -> lineno ∈ first_traj, run_df, view=true)
        push!(trajectories, trajectory)
        for i in 2:(length(traj_boundaries) - 1)
            trajectory = filter(:lineno => lineno -> lineno ∈ traj_boundaries[i]:(traj_boundaries[i+1] - 1), run_df, view=true)
            push!(trajectories, trajectory)
        end
        last_traj = last(traj_boundaries):only(describe(run_df, cols=:lineno, :max).max)
        trajectory = filter(:lineno => lineno -> lineno ∈ last_traj, run_df, view=true)
        push!(trajectories, trajectory)
    end
    return trajectories
end

function file_health_checks(output_df::DataFrame)
    if output_df[1, :].output != RUN_BEGIN
        @warn "File does not start with a run: " output_df[1, :].output
    end
end

function extract_run_metadata(runs)
    runs_metadata = DataFrame(finished_cleanly=Bool[], run_number=Int[], warnings=Any[], integrator=Any[], action=Any[], csw=Float64[], rng=Any[], β=Float64[], m0=Float64[], mu=Float64[], geometry=OffsetArray[])
    for run_number in 1:length(runs)
        metadata = Dict()
        metadata[:finished_cleanly] = false
        run = runs[run_number]

        if last(run).name == "SYSTEM" && last(run).output == "Process finalized."
            metadata[:finished_cleanly] = true
        end

        metadata[:run_number] = run_number
        metadata[:warnings] = filter(row -> row.name == "WARNING", run).output
        metadata[:integrator] = filter(row -> row.name == "INTEGRATOR", run).output
        metadata[:action] = filter(row -> row.name == "ACTION", run).output
        metadata[:csw] = parse(Float64, only(_extractor_only_one_matching_line(r"^Coefficient: reset to csw = (.*)$", "CLOVER", run, vital=true)))
        metadata[:rng] = filter(row -> row.name == "SETUP_RANDOM", run).output
        metadata[:β] = parse(Float64, _extractor_only_one_matching_line(r"^Monomial (.*): level = (.*), type = gauge, beta = (.*)$", "ACTION", run, vital=true)[3])
        tm_alt = _extractor_only_one_matching_line(r"^Monomial (.*): level = (.*), type = tm_alt, mass = (.*), mu = (.*), force_prec = (.*), mt_prec = (.*)$", "ACTION", run)
        hmc = _extractor_only_one_matching_line(r"^Monomial (.*): level = (.*), type = hmc, mass = (.*), force_prec = (.*), mt_prec = (.*)$", "ACTION", run)
        

        if(tm_alt != nothing)
            metadata[:m0] = parse(Float64, tm_alt[3])
            metadata[:mu] = parse(Float64, tm_alt[4])
        else
            metadata[:m0] = parse(Float64, hmc[3])
        end

        metadata[:geometry] = OffsetVector(parse.(Int, split(only(_extractor_only_one_matching_line(r"^Global size is (.*)$", "GEOMETRY_INIT", run, vital=true, multiple_but_unique_okay=true)), "x")), 0:3)
        
        
        push!(runs_metadata, [metadata[:finished_cleanly], metadata[:run_number], metadata[:warnings], metadata[:integrator], metadata[:action], metadata[:csw], metadata[:rng], metadata[:β], metadata[:m0], metadata[:mu], metadata[:geometry]]) 
    end
    return runs_metadata
end

function extract_wf_metadata(first_run)
    wf_metadata = Dict{Any, Any}()
    wf_metadata[:delta] = parse(Float64, only(_extractor_only_one_matching_line(r"^WF delta: (.*)$", "MAIN", first_run)))
    wf_metadata[:dt] = parse(Float64, only(_extractor_only_one_matching_line(r"^WF measurement interval dt : (.*)$", "MAIN", first_run)))
    
    return wf_metadata
end


function extract_trajectory_data(trajectories)
    TRAJ_MEASUREMENTS = [:accepted, :time, :mvm, :plaquette, :dS, :maxeig, :mineig, :adjoint_polyakov, :fundamental_polyakov]

    traj_data = DataFrame([Int[], Int[], Bool[], Union{Bool, Missing}[], Union{Float64, Missing}[], Union{Int64, Missing}[], Union{Float64, Missing}[], Union{Float64, Missing}[], Union{Float64, Missing}[], Union{Float64, Missing}[], Union{OffsetArray{Float64}, Missing}[], Union{OffsetArray{ComplexF64}, Missing}[], fill(Union{OffsetArray{Float64}, Missing}[], length(CORRELATORS))...], [:confno, :runno, :completed, TRAJ_MEASUREMENTS..., CORRELATORS...])
    for traj in trajectories
        data = Dict()
        data[:confno] = parse(Int, only(_extractor_only_one_matching_line(TRAJ_BEGIN_REGEX, "MAIN", traj, vital=true)))
        data[:runno] = only(union(traj.run_number))
        data[:completed] = true
        for measurement in cat(TRAJ_MEASUREMENTS, CORRELATORS, dims = 1)
            data[measurement] = missing
        end
        acc = _extractor_only_one_matching_line(ACCEPT_REJECT_REGEX, "HMC", traj)
        if isnothing(acc)
            data[:completed] = false
        else
            acc = only(acc)
            if acc == "accepted"
                data[:accepted] = true
            elseif acc == "rejected"
                data[:accepted] = false
            end
        end
        
        time = _extractor_only_one_matching_line(TRAJ_GENERATED_REGEX, "MAIN", traj)
	mvm = _extractor_only_one_matching_line(MVM_REGEX, "MAIN", traj)
        if isnothing(time)
            data[:completed] = false
	    data[:mvm] = missing
        else
            data[:time] = parse(Int, time[2]) + 1e-6*parse(Int, time[3])
            data[:mvm] = parse(Int, only(mvm))
        end

        plaquette = _extractor_only_one_matching_line(r"^Plaquette: (.*)$", "MAIN", traj)
        if isnothing(plaquette)
            data[:completed] = false
        else
            data[:plaquette] = parse(Float64, only(plaquette))
        end

        dS = _extractor_only_one_matching_line(DS_REGEX, "HMC", traj)
        if isnothing(dS)
            data[:completed] = false
        else
            data[:dS] = parse(Float64, dS[1])
        end

        maxeig = _extractor_only_one_matching_line(r"^Max_eig = (.*) \[MVM = .*\]$", "MaxH", traj)
        if !isnothing(maxeig)
            data[:maxeig] = parse(Float64, only(maxeig))
        end

        mineig = _extractor_only_one_matching_line(r"^Eig 0 = (.*)$", "LOWEIG", traj)
        if !isnothing(mineig)
            data[:mineig] = parse(Float64, only(mineig))
        end

        data[:adjoint_polyakov] = OffsetArray{Float64}(undef, 0:3)
        for i in 0:3
            adj_poly = _extractor_only_one_matching_line(r"^Polyakov direction " * string(i) * r" = (.*)$", "ADJ_POLYAKOV", traj)
            if isnothing(adj_poly)
                data[:adjoint_polyakov] = missing
                break
            else
                data[:adjoint_polyakov][i] = parse(Float64, only(adj_poly))
            end
        end
        
            
        data[:fundamental_polyakov] = OffsetArray{ComplexF64}(undef, 0:3)
        for i in 0:3
            fund_poly = _extractor_only_one_matching_line(r"^Polyakov direction " * string(i) * r" = (.*) (.*)$", "FUND_POLYAKOV", traj)
            if isnothing(fund_poly)
                data[:fundamental_polyakov] = missing
                break
            else
                data[:fundamental_polyakov][i] = parse(Float64, fund_poly[1]) + 1im * parse(Float64, fund_poly[2])
            end
        end

        corrs = []

        for corrname in CORRELATORS
            corrtemp = _extractor_only_one_matching_line(r"DEFAULT_SEMWALL TRIPLET " * string(corrname) * r"= (.*) $", "MAIN", traj)
            if isnothing(corrtemp)
                data[:completed] = false
                push!(corrs, missing)
            else
                c = parse.(Float64, split(only(corrtemp)))
                push!(corrs, OffsetVector(c, 0:(length(c) - 1)))
            end
        end

	push!(traj_data, [data[:confno], data[:runno], data[:completed], data[:accepted], data[:time], data[:mvm], data[:plaquette], data[:dS], data[:maxeig], data[:mineig], data[:adjoint_polyakov], data[:fundamental_polyakov], corrs...])
    end
    return traj_data
end

function extract_wf_trajectory_data(trajectories, metadata)
    obs = [:t, :E, :t2E, :Esym, :t2Esym, :TC, :W, :Wsym]
    traj_data = DataFrame([Int[], Union{Float64, Missing}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]], [:confno, :time, obs...])

    for traj in trajectories
        data = Dict()
        data[:confno] = parse(Int, only(_extractor_only_one_matching_line(WF_TRAJ_BEGIN_REGEX_RUN_NUMBER, "MAIN", traj, vital=true)))

        data[:time] = missing

        time = _extractor_only_one_matching_line(WF_TIMING_REGEX, "TIMING", traj)
        if !isnothing(time)
            data[:time] = parse(Int, time[1]) + 1e-6*parse(Int, time[2])
        end

        meas = []

        meas_lines = filter(:name => name -> name == "WILSONFLOW", traj).output
        for line in meas_lines
            push!(meas, parse.(Float64, split(only(match(WF_MEASUREMENT_REGEX, line).captures))))
        end
        meas = hcat(meas...)'
        data[:t] = meas[:, 1]
        data[:E] = meas[:, 2]
        data[:t2E] = meas[:, 3]
        data[:Esym] = meas[:, 4]
        data[:t2Esym] = meas[:, 5]
        data[:TC] = meas[:, 6]
        data[:W] = parent(d(data[:t2E], h = metadata[:dt])) .* data[:t][2 : end - 1]
        data[:Wsym] = parent(d(data[:t2Esym], h = metadata[:dt])) .* data[:t][2 : end - 1]

        push!(traj_data, [data[:confno], data[:time], data[:t], data[:E], data[:t2E], data[:Esym], data[:t2Esym], data[:TC], data[:W], data[:Wsym]])
    end
    
    return traj_data
end

function wilson_flow_same_tmax(trajectory_data::DataFrame)
    tmax = min(length.(trajectory_data.t)...)
    wmax = min(length.(trajectory_data.W)...)
    wf_meas = [:t, :E, :t2E, :Esym, :t2Esym, :TC]
    for i in 1:nrow(trajectory_data)
        for meas in wf_meas
            trajectory_data[i, meas] = trajectory_data[i, meas][1:tmax]
        end
        trajectory_data[i, :W] = trajectory_data[i, :W][1:wmax]
        trajectory_data[i, :Wsym] = trajectory_data[i, :Wsym][1:wmax]
    end
    return trajectory_data
end

function _extractor_only_one_matching_line(regex, name, df; vital=false, multiple_but_unique_okay=false)
    matching_line = filter([:name, :output] => (rowname, rowoutput) -> rowname == name && !isnothing(match(regex, rowoutput)), df).output
    if(multiple_but_unique_okay)
        matching_line = Set(matching_line)
    end
    try
        return match(regex, only(matching_line)).captures
    catch e
        if(vital)
            lineno = describe(df, cols=:lineno, :min).min
            @error "Error extracting data from dataframe starting at line number " lineno
            @show regex
            @show matching_line
            rethrow()
        end
    end
    return nothing
end

function _extractor_all_matching_lines(regex, name, df)
    return filter([:name, :output] => (rowname, rowoutput) -> rowname == name && !isnothing(match(regex, rowoutput)), df).output
end

function drop_runs_with_no_confs(data, run_metadata)
    return run_metadata[unique(data.runno), :]
end
    

function run_health_checks(run_metadata)

    #Check the geometry doesn't change
    if nrow((unique(select(run_metadata, :geometry)))) != 1
        @error "The geometry changed!"
        @show (unique(select(run_metadata, :geometry)))
    end


    
end

function find_nans(trajectory_data)

    #Check for NaNs
    nan_confs = []
    for corrname in CORRELATORS
        union!(nan_confs, filter(row -> !ismissing(row[corrname]) && any(isnan.(row[corrname])), trajectory_data).confno)
    end
    return nan_confs
    
end

function drop_missing_configurations(trajectories)
    return filter(row -> row.completed, trajectories)
end

function add_run_number_to_output_df(output_df::DataFrame)
    output_df[:, :run_number] .= 0
    run_boundaries = filter(:output => output -> output == RUN_BEGIN, output_df).lineno
    for i in 1:(length(run_boundaries) - 1)
        output_df[run_boundaries[i]:(run_boundaries[i+1] - 1), :run_number] .= i
    end
    output_df[last(run_boundaries):end, :run_number] .= length(run_boundaries)
    return output_df
end

function run_to_first_conf(data::DataFrame, run_number::Integer)
    confs_from_run = filter(:runno => runno -> runno == run_number, data)
    return only(describe(confs_from_run, cols=:confno, :min).min)
end

function post_process_correlators(trajectory_data)
    trajectory_data[:, :gk] = (trajectory_data[:, :g1] + trajectory_data[:, :g2] + trajectory_data[:, :g3])/3
    trajectory_data[:, :g5gk] = (trajectory_data[:, :g5g1] + trajectory_data[:, :g5g2] + trajectory_data[:, :g5g3])/3
    trajectory_data[:, :dg5_g0g5_re] = _corr_derivative(trajectory_data[:, :g5_g0g5_re])
    trajectory_data[:, :g5gk_folded] = _fold(trajectory_data[:, :g5gk])
    trajectory_data[:, :gk_folded] = _fold(trajectory_data[:, :gk])
    trajectory_data[:, :g5_folded] = _fold(trajectory_data[:, :g5])
    trajectory_data[:, :id_folded] = _fold(trajectory_data[:, :id])
    trajectory_data[:, :g5_g0g5_re_folded] = _fold(trajectory_data[:, :g5_g0g5_re], symm = false)
    trajectory_data[:, :dg5_g0g5_re_folded] = _corr_derivative(trajectory_data[:, :g5_g0g5_re_folded], last_point_antisymmetric = true)
    return trajectory_data
end

function extract_measurements_only(df::DataFrame; β, csw)
    masses = []
    meas_lines = _extractor_all_matching_lines(MEAS_REGEX, "MAIN", df)
    nconfs = parse(Int64, match(MEAS_REGEX, last(meas_lines)).captures[1])
    measurements = [Dict() for i in 1:nconfs]
    for line in meas_lines
        conf, m, corr, data = match(MEAS_REGEX, line).captures
        push!(masses, parse(Float64, m))
        data = parse.(Float64, split(data))
        measurements[parse(Int64, conf)][Symbol(corr)] = OffsetArray(data, 0:(length(data)-1))
    end

    mass = only(Set(masses))

    traj_data = DataFrame([Int[], fill(Union{OffsetArray{Float64}, Missing}[], length(CORRELATORS))...], [:confno, CORRELATORS...])

    nan_confs = find_nans(traj_data)

    for i in 1:nconfs
        push!(traj_data, [i; getindex.(Ref(measurements[i]), [CORRELATORS...])])
    end

    metadata = Dict()
    metadata[:geometry] = OffsetVector(parse.(Int, split(only(_extractor_only_one_matching_line(r"^Global size is (.*)$", "GEOMETRY_INIT", df, vital=true, multiple_but_unique_okay=true )), "x")), 0:3)
    metadata[:m0] = mass
    metadata[:β] = β
    metadata[:csw] = csw
    metadata[:nan_confs] = nan_confs
    metadata[:nconfs] = nrow(traj_data)

    traj_data = post_process_correlators(traj_data)
    
    return Ensemble(metadata, DataFrame([metadata]), traj_data, traj_data)
end


function load_ensemble(path::String; no_measurements = false)
    @debug "Loading the output file..."
    output_df = load_output_file_as_dataframe(path)
    @debug "Checking file health..."
    file_health_checks(output_df)
    @debug "Extracting runs..."
    output_df = add_run_number_to_output_df(output_df)
    runs = split_output_dataframe_into_runs(output_df, keep_runs_without_trajectories=true)
    @debug "Extracting run metadata..."
    run_metadata = extract_run_metadata(runs)
    @debug "Checking runs health..."
    run_health_checks(run_metadata)
    @debug "Extracting trajectories..."
    trajectories = split_run_dataframe_into_trajectories(runs)
    @debug "Extracting trajectory data..."
    trajectory_data = extract_trajectory_data(trajectories)
    @debug "Dropping missing configurations..."
    data = drop_missing_configurations(trajectory_data)
    run_metadata = drop_runs_with_no_confs(data, run_metadata)
    if !no_measurements
        @debug "Post-processing correlators"
        data = post_process_correlators(data)
    end
    @debug "Extracting global metadata..."
    global_metadata = extract_global_metadata(path, output_df, data, run_metadata)
    return Ensemble(global_metadata, run_metadata, data, data)
end

function load_wilsonflow(path::String)
    @debug "Loading the output file..."
    output_df = load_output_file_as_dataframe(path)
    @debug "Checking file health..."
    file_health_checks(output_df)
    @debug "Extracting runs..."
    output_df = add_run_number_to_output_df(output_df)
    runs = split_output_dataframe_into_runs(output_df, keep_runs_without_trajectories=true)
    @debug "Extracting metadata..."
    metadata = extract_wf_metadata(runs[1])
    @debug "Extracting trajectories..."
    trajectories = split_wf_dataframe_into_trajectories(runs)
    @debug "Extracting trajectory data..."
    trajectory_data = extract_wf_trajectory_data(trajectories, metadata)
    trajectory_data = wilson_flow_same_tmax(trajectory_data)
    return WilsonFlow(metadata, trajectory_data, trajectory_data)
end
    
function load_measurements(path::String; β, csw)
    @debug "Loading the output file..."
    output_df = load_output_file_as_dataframe(path)
    return extract_measurements_only(output_df, β = β, csw = csw)
end
