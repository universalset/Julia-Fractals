# Julia-Fractals
Julia language tools for generating fractal images and animations, mostly of Julia sets and related fractals.

This project is mostly for fun and amusement.  It may eventually evolve into being a proper package, but no promises.

## Getting started
See samples.jl or sample_notebook.ipynb for examples of use.

To set up:
Start terminal session, with root of the git repo as working directory.

Run the Julia REPL.
```
$ julia
```

Enter Pkg and activate.
```
julia> ] 
pkg> activate .
```
Now hit backspace to exit pkg.


To run the samples.jl script from here:
```
julia> include("samples.jl")
```

To start the notebook:
```
julia> using IJulia
julia> notebook(dir=pwd())
```

(Install Jupyter if needed.)

Open sample_notebook.ipynb from the interface.

Choose Cell->Run All

And then read, change, run, etc. cells as desired.