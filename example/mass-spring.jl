# Mass Spring

using ModelPredictiveControl

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

umin = [-0.5, -0.5, -0.5]
umax = [+0.5, +0.5, +0.5]
xmin = [-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0]
xmax = [+4.0, +4.0, +4.0, +4.0, +4.0, +4.0, +4.0, +4.0]

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
        stages.lb[i] .= umin
        stages.ub[i] .= umax
        stages.idx[i] .= 0:2
    elseif 1 < i < N+1
        stages.lb[i] .= [umin; xmin]
        stages.ub[i] .= [umax; xmax]
        stages.idx[i] .= 0:10
    else
        stages.lb[N+1] .= xmin
        stages.ub[N+1] .= xmax
        stages.idx[N+1] .= 3:10
    end

    # equality
    if i == 1
        stages.b[1] .= b .+ A*x₀
        stages.B[1] .= B
    end
    if 1 < i < N+1
        stages.A[i] .= A
        stages.b[i] .= b
        stages.B[i] .= B
    end

    # inequality
    if i == N+1
        stages.C[N+1] .= CN
    end
end

nbytes = getworkspacesize(HPMPCSolver, stages)
work = zeros(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work)
result = QuadraticProgramResult(stages)
solve!(result, stages, solver)
