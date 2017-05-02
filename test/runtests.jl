using ModelPredictiveControl
using Base.Test

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

@test length(stages.nx) == N+1
@test length(stages.nu) == N+1
@test length(stages.nb) == N+1
@test length(stages.nc) == N+1
@test stages.nx[1] == 0
@test all(stages.nx[2:end] .== nx)
@test stages.nu[N+1] == 0
@test all(stages.nu[1:end-1] .== nu)
@test stages.nb[1] == stages.nu[1]
@test stages.nb[end] == stages.nx[end]
@test all(stages.nb[2:end-1] .== nb)
@test stages.nc[end] == ncN
@test all(stages.nc[1:end-1] .== nc)

@test size(stages.A[1]) == (nx, 0)
@test all(map(size, stages.A[2:end]) .== (nx, nx))
@test all(map(size, stages.B) .== (nx, nu))
@test all(map(size, stages.b) .== (nx,))
@test size(stages.Q[1]) == (0, 0)
@test all(map(size, stages.Q[2:end]) .== (nx, nx))
@test size(stages.S[1]) == (nu, 0)
@test all(map(size, stages.S[2:end]) .== (nu, nx))
@test all(map(size, stages.R) .== (nu, nu))
@test size(stages.q[1]) == (0,)
@test all(map(size, stages.q[2:end]) .== (nx,))
@test all(map(size, stages.r) .== (nu,))
@test size(stages.lb[1]) == (nu,)
@test size(stages.ub[1]) == (nu,)
@test size(stages.lb[end]) == (nx,)
@test size(stages.ub[end]) == (nx,)
@test all(map(size, stages.lb[2:end-1]) .== (nb,))
@test all(map(size, stages.ub[2:end-1]) .== (nb,))
@test size(stages.idx[1]) == (nu,)
@test size(stages.idx[end]) == (nx,)
@test all(map(size, stages.idx[2:end-1]) .== (nb,))
@test size(stages.lc[end]) == (ncN,)
@test size(stages.uc[end]) == (ncN,)
@test all(map(size, stages.lc[1:end-1]) .== (nc,))
@test all(map(size, stages.uc[1:end-1]) .== (nc,))
@test size(stages.C[1]) == (nc, 0)
@test size(stages.C[end]) == (ncN, nx)
@test all(map(size, stages.C[2:end-1]) .== (nc, nx))
@test all(map(size, stages.D) .== (nc, nu))

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
    @test status == 0
    #TODO add more tests
    U[:, k] .= result.u[1]
    X[:, k+1] = A*X[:, k] + B*U[:, k]
end
