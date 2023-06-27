#!julia

confdir = ARGS[1]

if first(confdir) != '/'
	@error "Use absolute paths only."
	exit(1)
end

if last(confdir) != '/'
	confdir *= '/'
end

confs = readdir(ARGS[1])
n = []
lines = []
for conf in confs
	c = match(r"^(.*n)([0-9]+)$", conf).captures
	push!(n, parse(Int64, c[2]))
	push!(lines, c[1])
end

line = only(Set(lines))
sort!(n)

open(ARGS[2], "w") do f
	for i in n
		write(f, confdir*line*string(i)*"\n")
	end
end
