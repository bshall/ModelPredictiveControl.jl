using ModelPredictiveControl
using ControlSystems
using JLD

N = 20
nx = 6
nu = 4
nb = 4
nc = 0
ncN = 0

data = load("example/suspension.jld")
A = data["A"]
B = data["B"]
Bw = data["Bw"]
Cd = data["Cd"]
Dd = data["Dd"]

Q = 50*eye(nx)
R = eye(nu)
P = dare(A, B, Q, R)

umin = [-0.04, -0.04, -0.04, -0.04]
umax = [+0.04, +0.04, +0.04, +0.04]

# road
n = 3600
t = 0:0.005:0.005n-0.005
road = zeros(length(t))
road[161:180] = 0.005*(1:20)
road[181:200] = 0.1 - 0.005*(1:20)

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
        stages.idx[i] .= 0:nu-1
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
solver = HPMPCSolver(work, 1e-8, 20, 50.0, false)
result = QuadraticProgramResult(stages)

steps = 720
X = zeros(nx, steps+1)
X[:, 1] = zeros(nx)
Y = zeros(3, steps)
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
