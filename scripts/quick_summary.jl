#!julia
include(Base.source_dir()*"/../src/parser.jl")
ens = load_ensemble(ARGS[1])

if !isempty(ens.global_metadata[:integrator_changes])
    data_since_last_change = filter(:confno => confno -> confno >= max(ens.global_metadata[:integrator_changes]...), ens.data)
else
    data_since_last_change = ens.data
end
   

@info "Latest configuration: "
println(ens.data[end, [:confno, :accepted, :time, :plaquette]])
println()

@info "Statistics since last integrator change: "
println("Number of confs since last change: ", nrow(data_since_last_change))
println(describe(data_since_last_change, :mean, :median, cols=[:accepted, :time, :plaquette]))
println()
        
@info "Global metadata:"
@show ens.global_metadata
