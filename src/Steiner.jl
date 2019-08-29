module Steiner

using DelimitedFiles, LinearAlgebra, StaticArrays

import HomotopyContinuation, DynamicPolynomials
const DP = DynamicPolynomials
const HC = HomotopyContinuation

include("solve.jl")
include("pipe.jl")
include("server.jl")

const startconics_startsolutions = begin
    startconics = SVector{30}(readdlm(joinpath(@__DIR__, "..", "deps", "startconics.txt"), '\t', ComplexF64))
    allsolutions_real = readdlm(joinpath(@__DIR__, "..", "deps", "allsolutions_real.txt"), '\t', Float64, '\n')
    allsolutions_imag = readdlm(joinpath(@__DIR__, "..", "deps", "allsolutions_imag.txt"), '\t', Float64, '\n')
    allsolutions = [[complex(allsolutions_real[i, j], allsolutions_imag[i,j]) for j=1:15] for i=1:3264]
    startconics, allsolutions
end
const trackers = [assemble_tracker()]


function __init__()
    for i in 2:Threads.nthreads()
        push!(trackers, deepcopy(trackers[1]))
    end
end

end # module
