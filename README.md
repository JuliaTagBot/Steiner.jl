# Steiner

Compute all conics tangents to five given tonics.

## Installation

Start a Julia session and enter the `pkg`-mode by pressing `]`:
```julia
pkg> add https://github.com/JuliaHomotopyContinuation/Steiner.jl.git
```
alternatively, you can do
```julia
julia> using Pkg; Pkg.add("https://github.com/JuliaHomotopyContinuation/Steiner.jl.git");
```

We provide a web interface to compute the conics tangent. You can start this with
```julia
julia> using Steiner;
julia> start_server();
```

By default you then should be able to go to [http://localhost:3264](http://localhost:3264)
and compute your conics.
