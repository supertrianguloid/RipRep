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

mutable struct Ensemble
    global_metadata::DataFrame
    run_metadata::DataFrame
    data::DataFrame
    analysis::DataFrame
end

mutable struct WilsonFlow
    metadata::DataFrame
    data::DataFrame
end

ensembledir() = datadir() * "/ensembles"
reportsdir() = projectdir() * "/reports"

output_file_location(ens_no) = ensembledir() * "/" * GLOBAL_SIMS[ens_no, :path] * "/out_0"
runs_ignore_file(ens_no) = ensembledir() * "/" * GLOBAL_SIMS[ens_no, :path] * "/RipRep/ignore_runs"

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

function load_ensemble(ensemble_no::Integer)
    g_meta, r_meta, data = parse_output_file(ensemble_no)

    g_meta = hcat(DataFrame(GLOBAL_SIMS[ensemble_no, :]), g_meta)

    e = Ensemble(g_meta, r_meta, data, data)
    e.global_metadata.IntegratorChanges = [_check_integrator(e)]
    e.global_metadata.Acceptance = [_measure_acceptance(e)]
    return e
end

function load_ensemble(path::String)
    g_meta, r_meta, data = parse_output_file(ensemble_no)

    g_meta = hcat(DataFrame(GLOBAL_SIMS[ensemble_no, :]), g_meta)

    e = Ensemble(g_meta, r_meta, data, data)
    e.global_metadata.IntegratorChanges = [_check_integrator(e)]
    e.global_metadata.Acceptance = [_measure_acceptance(e)]
    return e
end

function summarise_ensembles(method = :show)
    df = DataFrame()
    for i in 1:nrow(GLOBAL_SIMS)
        println("Parsing ensemble $i...")
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