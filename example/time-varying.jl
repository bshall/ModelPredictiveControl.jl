using ModelPredictiveControl

N = 15
nx = 2
nu = 1
nb = 1
nc = 0
ncN = 0

A₁ = [0.7115 -0.6; 0.6 0.8853]
B₁ = [0.02713, 0.0573]

A₂ = [0.9 0.5; 0.5 1.0]
B₂ = [0, 0.6667]

A₃ = [0.7115 -0.5; 0.5 1.0]
B₃ = [0.5, 0.01]

A₄ = [0.0 0.9; -1.0 0.0]
B₄ = [0.0 0.2]

Q = 10*eye(nx)
R = eye(nu)

u₁min = [-3.0]
u₁max = [+5.0]

u₂min = [-5.5]
u₂max = [+5.5]

u₃min = [-3.0]
u₃max = [+5.0]

u₄min = [-0.45]
u₄max = [+4.5]

s₁, s₂, s₃, s₄ = 4, 8, 3, 5

stages = MultistageProblem(N, nx, nu, nb, nc, ncN)

for i = 1:N+1
    # costs
    if i == 1
        stages.R[1] .== R
    elseif 1 < i < N+1
        stages.Q[i] .== Q
        stages.R[i] .== R
    else
        stages.Q[N+1] .== Q
    end
end

nbytes = getworkspacesize(HPMPCSolver, stages)
work = zeros(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work)
result = QuadraticProgramResult(stages)

steps = 40
X = zeros(nx, steps+1)
X[:, 1] = [-4.0, 2.0]
U = zeros(nu, steps)
for k = 1:steps
    stages.b[1] .==
end
