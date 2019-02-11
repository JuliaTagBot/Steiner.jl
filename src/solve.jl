using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, StaticArrays
import ProjectiveVectors
using DelimitedFiles

function steiner_system()
    @polyvar x[1:2] a[1:5] c[1:6] y[1:2, 1:5]

    #tangential conics
    f = a[1]*x[1]^2 + a[2]*x[1]*x[2] + a[3]*x[2]^2 + a[4]*x[1]  + a[5]*x[2] + 1;
    ∇ = differentiate(f, x)
    #5 conics
    g = c[1]*x[1]^2 + c[2]*x[1]*x[2] + c[3]*x[2]^2 + c[4]*x[1]  + c[5]*x[2] + c[6];
    ∇_2 = differentiate(g, x)
    #the general system
    #f_a_0 is tangent to g_b₀ at x₀
    function Incidence(f,a₀,g,b₀,x₀)
        fᵢ = f(x=>x₀, a=>a₀)
        ∇ᵢ = [∇ᵢ(x=>x₀, a=>a₀) for ∇ᵢ in ∇]
        Cᵢ = g(x=>x₀, c=>b₀)
        ∇_Cᵢ = [∇ⱼ(x=>x₀, c=>b₀) for ∇ⱼ in ∇_2]

        [fᵢ; Cᵢ;  det([∇ᵢ ∇_Cᵢ])]
    end

    #Track the solutions elsewhere
    @polyvar v[1:6, 1:5]
    F = vcat(map(i -> Incidence(f,a,g, v[:,i], y[:,i]), 1:5)...)
    F, vec(v)
end

const startconics_startsolutions = begin
    startconics = SVector{30}(readdlm(joinpath(@__DIR__, "..", "deps", "startconics.txt"), '\t', ComplexF64))
    allsolutions_real = readdlm(joinpath(@__DIR__, "..", "deps", "allsolutions_real.txt"), '\t', Float64, '\n')
    allsolutions_imag = readdlm(joinpath(@__DIR__, "..", "deps", "allsolutions_imag.txt"), '\t', Float64, '\n')
    allsolutions = [[complex(allsolutions_real[i, j], allsolutions_imag[i,j]) for j=1:15] for i=1:3264]
    startconics, allsolutions
end

function assemble_tracker()
    startconics, startsolutions = startconics_startsolutions
    F, parameters = steiner_system()
    pathtracker(F, startsolutions, parameters=parameters, p₁=startconics, p₀=(@SVector randn(30)),
        tol=1e-6,
        refinement_tol=1e-10)
end

const tracker = assemble_tracker()

function setup_parameters!(homotopy, p₁, p₀)
    H = basehomotopy(homotopy)
    if !(H isa ParameterHomotopy)
       error("Base homotopy is not a ParameterHomotopy")
    end
    set_parameters!(H, (p₁, p₀))
end


"""
    solve_conics(M::Matrix)

The conics are described as 5 × 5 matrix.
Each column describes a normalized conic, i.e., with constant term 1.

Output is a vector of normalized conics.
"""
function solve_conics(M::Matrix)
    tangential_conics = Vector{Float64}[]
    tangential_points = Vector{Float64}[] # Stores [x1,y1,x2,y2,...,x5,y5]
    complex_solutions = Vector{ComplexF64}[]
    tangential_indices = Int[]
    compute_time = @elapsed solve_conics(M) do x, k
        if HomotopyContinuation.isrealvector(x, 1e-8)
            push!(tangential_conics, real.(x[1:5]))
            push!(tangential_points, real.(x[6:15]))
            push!(tangential_indices, k)
        end
        push!(complex_solutions, x[1:5])
    end
    nreal = length(tangential_conics)

    Dict("tangential_conics" => tangential_conics,
         "tangential_points" => tangential_points,
         "tangential_indices" => tangential_indices,
         "nreal" => nreal,
         "complex_solutions" => Dict("real" => real.(complex_solutions),
                                     "imag" => imag.(complex_solutions)))
end

function solve_conics(output!::F, M::Matrix) where {F<:Function}
    p₁, startsolutions = startconics_startsolutions
    p₀ = SVector{30}(complex.(M))
    setup_parameters!(tracker.homotopy, p₁, p₀)
    solve_conics(output!, tracker, startsolutions)
    nothing
end

function solve_conics(output!::F, pathtracker::PathTracker, starts) where {F<:Function}
    for (k, s) in enumerate(starts)
        ret = track!(pathtracker, s, 1.0, 0.0)
        if ret == PathTrackerStatus.success
            x = ProjectiveVectors.affine_chart(currx(pathtracker))
            output!(x, k)
        end
    end
    nothing
end


function handle_solve_conics(conics_input)
    M = ones(6, 5)
    @info "Compute conics: $conics_input"
    for j in 1:5, i in 1:5
        M[i, j] = conics_input[j][i]
    end
    @info "M: $M"
    solve_conics(M)
end
