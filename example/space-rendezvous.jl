using ModelPredictiveControl

# Scenario constants
const G = 6.673e-11           # Gravitational constant
const Mmars = 6.4191e23       # Mass of central body
const ecc = 0.204410          # Eccentricity
const a = 4643e3              # Semimajor axis a

# Calculations
const μ = G*Mmars                         # Gravitational parameters
const n = sqrt(μ/a^3)                     # Mean anomaly rate
const p = a*(1 - ecc^2)                   # Semi-latus rectum
const rp = p/(1 + ecc)                    # Periapsis
const v₀ = sqrt(2μ)*sqrt(1/rp - 1/(2*a))  # Velocity at periapsis
const h = rp*v₀                           # Angular momentum
const ya = h/(p^2)                        # Useful constant
const a₁ = sqrt(1-ecc)
const a₂ = sqrt(1+ecc)
const a₃ = ecc*sqrt(1-ecc^2)

function propagate_anomaly(M₀, Δt)
    M = M₀ + n*Δt
    E = M
    for k = 1:6
        ɛ = E - ecc*sin(E) - M
        δ = 1 - ecc*cos(E)
        E = E - ɛ/δ
    end
    ν = atan2(sin(E)*sqrt(1-ecc^2), cos(E)-ecc)
    return ν
end

function yamanaka_ankersen_model(ν, Δt)
    ρ = 1 + ecc*cos(ν)
    ϕinv = [1-ecc^2  0        3*ecc*ρ*sin(ν)*(1/ρ + 1/ρ^2)  -ecc*ρ*sin(ν)*(1+1/ρ)  0        -ecc*ρ*cos(ν)+2;
            0        1-ecc^2  0                             0                      0        0;
            0        0        -3*ρ*sin(ν)*(1/ρ+ecc^2/ρ^2)   ρ*sin(ν)*(1+1/ρ)       0        ρ*cos(ν)-2*ecc;
            0        0        -3*(cos(ν) + ecc)             ρ*cos(ν)*(1+1/ρ)+ecc   0        -ρ*sin(ν);
            0        0        0                             0                      1-ecc^2  0;
            0        0        3*ρ+ecc^2-1                   -ρ^2                   0        ecc*ρ*sin(ν)]/(1-ecc^2)

    M = 2*atan2(a₁*tan(ν/2), a₂) - (a₃*sin(ν))/(1+ecc*cos(ν))
    νnew = propagate_anomaly(M, Δt)

    ρ = 1 + ecc*cos(νnew)
    ds = cos(νnew)+ecc*cos(2*νnew)
    dc = -(sin(νnew)+ecc*sin(2νnew))
    ϕ = [1  0             -ρ*cos(νnew)*(1+1/ρ)  ρ*sin(νnew)*(1+1/ρ)  0             3*ρ^2*ya*Δt;
         0  cos(νnew-ν)   0                     0                    sin(νnew-nu)  0;
         0  0             ρ*sin(νnew)           ρ*cos(νnew)          0             2-3*ecc*ρ*sin(νnew)*ya*Δt;
         0  0             2*ρ*sin(νnew)         2*ρ*cos(νnew)-ecc    0             3*(1-2*ecc*ρ*sin(νnew)*ya*Δt);
         0  -sin(νnew-ν)  0                     0                    cos(νnew-ν)   0;
         0  0             ds                    dc                   0             -3*ecc*(ds*ya*Δt + sin(νnew)/ρ)]

    transk = [eye(3) zeros(3, 3);
              eye(3) eye(3)]
    itransk = [eye(3) zeros(3, 3);
               eye(3) eye(3)]
    return itransk*ϕnew*ϕinv*transk
end

N = 30
Δt = 200
nx = 6
nu = 3
nb = 3
nc = 3

Q = diagm([0.001, 0.001, 0.001, 0.000, 0.000, 0.000])
R = diagm([15000, 15000, 15000])

umin = [-5.0, -5.0, -5.0]
umax = [+5.0, +5.0, +5.0]

coneangle = 20*pi/180
Hx = [  1.0     0.0  0.0  0.0  0.0  0.0;
      coneangle 0.0  1.0  0.0  0.0  0.0;
      coneangle 0.0 -1.0  0.0  0.0  0.0]
hx = [0.0, 0.0, 0.0]

xs = [-1000, 0.0, 0.0, 0.0, 0.0, 0.0]
us = [0.0, 0.0, 0.0]

stages = MultistageProblem(N, nx, nu, )

for i = 1:N+1
    # costs
    if i == 1
        stages.R[1] .= R
        stages.r[1] .= -R*us
    elseif 1 < i < N+1
        stages.Q[i] .= Q
        stages.R[i] .= R
        stages.q[i] .= -Q*xs
        stages.r[i] .= -R*us
    else
        stages.Q[N+1] .= Q
        stages.q[N+1] .= -Q*xs
    end

    # bounds
    if i < N+1
        stages.lb[i] .= umin
        stages.ub[i] .= umax
        stages.idx[i] .= 0:nu-1
    end
end

nbytes = getworkspacesize(HPMPCSolver, stages)
work = zeros(cld(nbytes, sizeof(Float64)))
solver = HPMPCSolver(work, 1e-8, 40, 10.0, false)
result = QuadraticProgramResult(stages)

steps = 6000
X = zeros(nx, steps+1)
X[:, 1] = [-15e3, 0.0, 0.0, 0.0, 0.0, 0.0]
U = zeros(nu, steps)
ν = 0
for k = 1:steps
    if mod(k-1, Δt) == 0
        νcopy = ν
        A, νcopy = dynamics(νcopy, )
        B = A*[zeros(3); eye(3)]
        stages.b[1] .= A*X[:, k]
        stages.B[1] .= B
        for i = 2:N
            A, νcopy = dynamics(νcopy, )
            B = A*[zeros(3); eye(3)]
            stages.A[i] .= A
            stages.B[i] .= B
        end
        status = solve!
        if status == 0
            U[:, k] .= result.u[1]
        else
            error("Some problem in solver")
        end
    end
    A, ν = dynamics(ν, )
    X[:, k+1] = A*X[:, k] + B*U[:, k]
end
