using ModelPredictiveControl
using ControlSystems

N = 10
nx = 2
nu = 1
nb = 1
nc = 0
ncN = 0

A = [0.7115 -0.4345; 0.4345 0.8853]
B = [1.0, 1.0]
Bw = [1.0, 1.0]

Q = 10*eye(nx)
R = eye(nu)
P = dare(A, B, Q, R)

umin = [-1.8]
umax = [+1.8]

# generate disturbance
n = 100
road = zeros(n)
road[16:25] = 0.2*(1:10)
road[26:35] = 2 - 0.2*(1:10)

stages = MultistageProblem(N, nx, nu, nb, nc, ncN)

for i = 1:N+1
    # costs
    if i == 1
        stages.R[1] .= R
    elseif 1 < i < N+1
        stages.Q[i] .= Q
        stages.R[i] .= R
    else
        stages.Q[N+1] .= P
    end

    # bounds
    if i < N+1
        stages.lb[i] .= umin
        stages.ub[i] .= umax
    end

    # equality
    if i == 1
        stages.B[1] .= B
    end
    if 1 < i < N+1
        stages.A[i] .= A
        stages.B[i] .= B
    end
end

nbytes = getworkspacesize(HPMPCSolver, stages)
work = zeros(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work, 1e-8, 20, 10.0, false)
result = QuadraticProgramResult(stages)

steps = 60
X = zeros(nx, steps+1)
X[:, 1] = [0.0, 0.0]
U = zeros(nu, steps)
for k = 1:steps
    stages.b[1] .= Bw.*road[k] .+ A*X[:, k]
    for i = 2:N
        stages.b[i] .= Bw.*road[k+i-1]
    end
    status = solve!(result, stages, solver)
    if status == 0
        U[:, k] .= result.u[1]
    else
        error("Some problem in solver")
    end
    X[:, k+1] = A*X[:, k] + [B Bw]*[U[:, k]; road[k]]
end
