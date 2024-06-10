#!julia

JOBFILE = "job"
EXEFILE = "compute_loops"
CONFFILE = ARGS[1]

for conf in eachline(CONFFILE)
	confno = lpad(last(split(conf, 'n')), 4, "0")
	mkdir(confno)
	cp(JOBFILE, "$confno/job")
	cp(EXEFILE, "$confno/compute_loops")
	cd(confno)
	write("list", conf)
	run(`sbatch job`)
	cd("../")
end

