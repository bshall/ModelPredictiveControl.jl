module ModelPredictiveControl

export QuadraticProgram,
       QuadraticProgramResult,
       HPMPCSolver,
       getworkspacesize,
       solve!

depsjl = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(depsjl)
    include(depsjl)
else
    error("ModelPredictiveControl is not properly installed.
           Please run Pkg.build(\"ModelPredictiveControl\") and restart julia")
end

abstract AbstractSolver

include("QuadraticProgram.jl")
include("HPMPC.jl")

end # module
