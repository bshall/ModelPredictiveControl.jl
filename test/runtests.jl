using ModelPredictiveControl
using Base.Test

# Example 1

N = 50
nx = 2
nu = 1
nb = 1
nc = 0
ncN = 0

A = [1.0 0.1; 0.0 1.0]
B = [0.0, 0.1]
b = [0.0, 0.0]
Q = diagm([10.0, 1.0])
lb = [-1.0]
ub = [+1.0]

x₀ = [-1.0, 0.0]

stages = MultistageProblem(N, nx, nu, nb, nc, ncN)

for i = 1:N+1
    # costs
    if 1 < i <= N+1
        stages.Q[i] .= Q
    end

    # bounds
    if 1 <= i < N+1
        stages.lb[i] .= lb
        stages.ub[i] .= ub
    end

    # equality
    if i == 1
        stages.b[1] .= b .+ A*x₀
    elseif 1 < i < N+1
        stages.A[i] .= A
        stages.b[i] .= b
        stages.B[i] .= B
    end
end

nbytes = getworkspacesize(HPMPCSolver, stages)
work = Vector{Float64}(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work)
result = QuadraticProgramResult(stages)
solve!(result, stages, solver)

# Example 2

N = 15
nx = 2
nu = 1


# Mass Spring

N = 15
nx = 8
nu = 3
nb = 11
nc = 0
ncN = 4

Q = eye(nx)
R = 2*eye(nu)
S = zeros(nu, nx)
q = 0.1*ones(nx)
r = 0.2*ones(nu)

Δt = 0.5
Ac = [zeros(4, 4) eye(4); diagm(-2*ones(4)) + diagm(ones(3), -1) + diagm(ones(3), 1) zeros(4, 4)]
Bc = [zeros(4, 3); eye(3); zeros(1, 3)]
M = expm([Δt*Ac Δt*Bc; zeros(3, 11)])
A = M[1:8, 1:8]
B = M[1:8, 9:end]
b = 0.1*ones(nx)

lb = [-0.5, -0.5, -0.5, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0]
ub = [+0.5, +0.5, +0.5, +4.0, +4.0, +4.0, +4.0, +4.0, +4.0, +4.0, +4.0]

CN = [eye(4) zeros(4, 4)]

x₀ = [2.5, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

stages = MultistageProblem(N, nx, nu, nb, nc, ncN)

for i=1:N+1
    # costs
    if i == 1
        stages.R[1] .= R
        stages.r[1] .= r .+ S*x₀
    elseif 1 < i < N+1
        stages.Q[i] .= Q
        stages.R[i] .= R
        stages.S[i] .= S
        stages.q[i] .= q
        stages.r[i] .= r
    else
        stages.Q[N+1] .= Q
        stages.q[N+1] .= q
    end

    # bounds
    if i == 1
        stages.lb[i] .= lb[1:3]
        stages.ub[i] .= ub[1:3]
        stages.idx[i] .= 0:2
    elseif 1 < i < N+1
        stages.lb[i] .= lb
        stages.ub[i] .= ub
        stages.idx[i] .= 0:10
    else
        stages.lb[N+1] .= lb[4:end]
        stages.ub[N+1] .= ub[4:end]
        stages.idx[N+1] .= 3:10
    end

    # equality
    if i == 1
        stages.b[1] .= b .+ A*x₀
    elseif 1 < i < N+1
        stages.A[i] .= A
        stages.b[i] .= b
        stages.B[i] .= B
    end

    # inequality
    if i == N+1
        stages.C[N+1] .= CN
    end
    # if i == 1
    #     stages.lc[1] .= lc .- C*x₀
    #     stages.uc[1] .= uc .- C*x₀
    # elseif 1 < i < N+1
    #     stages.lc[i] .= lc
    #     stages.uc[i] .= uc
    # else
    #     stages.l
    #     stages.C[i] .= C
    # end
end

nbytes = getworkspacesize(HPMPCSolver, stages)
work = Vector{Float64}(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work)
result = QuadraticProgramResult(stages)
solve!(result, stages, solver)

result.u
result.x
solver.res
