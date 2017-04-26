using ModelPredictiveControl
using Base.Test

N = 5
nx = 2
nu = 1
nb = 0
nc = 0

A = [0.0 1.0; 0.0 0.0]
B = [1.0; 1.0]
Q = diagm([1.0, 1.0])
R = diagm([1.0])

x₀ = [1.0, 1.0]

qp = QuadraticProgram(N, nx, nu, nb, nc)
qp.lb[1] .= x₀
qp.ub[1] .= x₀
for i = 1:qp.N
    qp.A[i] .= A
    qp.B[i] .= B
    qp.Q[i] .= Q
    qp.R[i] .= R
end
qp.Q[qp.N+1] .= Q

nbytes = getworkspacesize(HPMPCSolver, qp)
work = Vector{Float64}(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work)
result = QuadraticProgramResult(qp)
