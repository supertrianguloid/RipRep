# RipRep

Reproducible *HiRep* analysis, including a reasonably sophisticated parser.

To reproduce this project, do the following:

0. Download this code base. Put raw *HiRep* files into `data/`
1. Open a Julia console and run:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
