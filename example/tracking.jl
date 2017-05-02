# Reference Tracking

using ModelPredictiveControl

N = 10
nu = 2
nx = 2
nb = 4
nc = 0
ncN = 0

Q = eye(nx)
R = eye(nu)

umin = [-0.8, -0.8]
umax = [+1.3, +1.3]
xmin = [-5.0, -5.0]
xmax = [+5.0, +5.0]

A = [1.1 1.0; 0.0 1.0]
B = [0.8 0.1; 0.3 0.8]

stages = MultistageProblem(N, nx, nu, nb, nc, ncN)

for i = 1:N+1
    # costs
    if i == 1
        stages.R[1] .= R
    elseif 1 < i < N+1
        stages.Q[i] .= Q
        stages.R[i] .= R
    else
        stages.Q[N+1] .= Q
    end

    # bounds
    if i == 1
        stages.lb[1] .= umin
        stages.ub[1] .= umax
        stages.idx[1] .= 0:nu-1
    elseif 1 < i < N+1
        stages.lb[i] .= [umin; xmin]
        stages.ub[i] .= [umax; xmax]
        stages.idx[i] .= 0:nu+nx-1
    else
        stages.lb[N+1] .= xmin
        stages.ub[N+1] .= xmax
        stages.idx[N+1] .= nu:nu+nx-1
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
solver = HPMPCSolver(work)
result = QuadraticProgramResult(stages)

steps = 35
X = zeros(nx, steps+1)
X[:, 1] = [-5.0, 4.0]
Xref = zeros(nx, steps+1)
Xref[1, 14:end] = 2
Xref[2, 22:end] = -1
U = zeros(nu, steps)
Uref = zeros(nu, steps)
for k = 1:steps
    Uref[:, k] = inv(B)*(eye(nx) - A)*Xref[:, k]
    stages.b[1] .= A*X[:, k]
    for i = 1:N+1
        if i == 1
            stages.r[1] .= -R'*Uref[:, k]
        elseif 1 < i < N+1
            stages.r[i] .= -R'*Uref[:, k]
            stages.q[i] .= -Q'*Xref[:, k]
        else
            stages.q[N+1] .= -Q'Xref[:, k]
        end
    end
    status = solve!(result, stages, solver)
    if status == 0
        U[:, k] .= result.u[1]
    else
        error("Some problem in solver")
    end
    X[:, k+1] = A*X[:, k] + B*U[:, k]
end
