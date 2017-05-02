using ModelPredictiveControl

# Scenario constants
const G = 6.673e-11           # Gravitational constant
const Mmars = 6.4191e23       # Mass of central body
const ecc = 0.204410          # Eccentricity
const a = 4643e3              # Semimajor axis a

# Calculations
const μ = G*Mmars                         # Gravitational parameters
const n = sqrt(μ/a^3)                     # Mean anomaly rate
const p = a*(1 - ecc^2)                    # Semi-latus rectum
const rp = p/(1 + ecc)                    # Periapsis
const v₀ = sqrt(2μ)*sqrt(1/rp - 1/(2*a))  # Velocity at periapsis
const h = rp*v₀                           # Angular momentum
const yaₖ = h/(p^2)                       # Useful constant

X1 = [sqrt(1 - ecc), sqrt(1 + ecc)]
X2 = ecc*sqrt(1 - ecc^2)

function kepler_mean_anomaly_rate(μ, a)
    
end

function dynamics(ν₀, Δt)
    ρ = 1 + ecc*cos(ν₀)

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
for k = 1:steps
    if mod(k-1, Δt) == 0

    end
end
