type MultistageProblem
    N::Cint
    nx::Vector{Cint}
    nu::Vector{Cint}
    nb::Vector{Cint}
    nc::Vector{Cint}

    A::Vector{Matrix{Cdouble}}
    B::Vector{Matrix{Cdouble}}
    b::Vector{Vector{Cdouble}}

    Q::Vector{Matrix{Cdouble}}
    S::Vector{Matrix{Cdouble}}
    R::Vector{Matrix{Cdouble}}
    q::Vector{Vector{Cdouble}}
    r::Vector{Vector{Cdouble}}

    lb::Vector{Vector{Cdouble}}
    ub::Vector{Vector{Cdouble}}
    idx::Vector{Vector{Cint}}

    lc::Vector{Vector{Cdouble}}
    uc::Vector{Vector{Cdouble}}
    C::Vector{Matrix{Cdouble}}
    D::Vector{Matrix{Cdouble}}

    MultistageProblem(N::Int) = new(N)
end

function MultistageProblem(N::Int, nx::Int, nu::Int, nb::Int, nc::Int, ncN::Int)
    nxvect = fill(nx, N+1)
    nuvect = fill(nu, N+1)
    nbvect = fill(nb, N+1)
    ncvect = fill(nc, N+1)
    nxvect[1] = 0
    nuvect[N+1] = 0
    nbvect[1] = (nb < nu) ? nb : nu
    nbvect[N+1] = (nb > nu) ? nb - nu : 0
    ncvect[N+1] = ncN

    stages = MultistageProblem(N, nxvect, nuvect, nbvect, ncvect)
    return stages
end

function MultistageProblem(N::Int, nx::Vector{Int}, nu::Vector{Int},
                           nb::Vector{Int}, nc::Vector{Int})
    length(nx) == N+1 || error()
    length(nu) == N+1 || error()
    length(nb) == N+1 || error()
    length(nc) == N+1 || error()

    stages = MultistageProblem(N)

    stages.nx = nx
    stages.nu = nu
    stages.nb = nb
    stages.nc = nc

    stages.A = [zeros(Cdouble, nx[i+1], nx[i]) for i=1:N]
    stages.B = [zeros(Cdouble, nx[i+1], nu[i]) for i=1:N]
    stages.b = [zeros(Cdouble, nx[i+1]) for i=1:N]

    stages.Q = [zeros(Cdouble, nx[i], nx[i]) for i=1:N+1]
    stages.S = [zeros(Cdouble, nu[i], nx[i]) for i=1:N]
    stages.R = [zeros(Cdouble, nu[i], nu[i]) for i=1:N]
    stages.q = [zeros(Cdouble, nx[i]) for i=1:N+1]
    stages.r = [zeros(Cdouble, nu[i]) for i=1:N]

    stages.lb = [zeros(Cdouble, nb[i]) for i=1:N+1]
    stages.ub = [zeros(Cdouble, nb[i]) for i=1:N+1]
    stages.idx = [zeros(Cint, nb[i]) for i=1:N+1]

    stages.lc = [zeros(Cdouble, nc[i]) for i=1:N+1]
    stages.uc = [zeros(Cdouble, nc[i]) for i=1:N+1]
    stages.C = [zeros(Cdouble, nc[i], nx[i]) for i=1:N+1]
    stages.D = [zeros(Cdouble, nc[i], nu[i]) for i=1:N]

    return stages
end

type QuadraticProgramResult
    x::Vector{Vector{Cdouble}}
    u::Vector{Vector{Cdouble}}
    π::Vector{Vector{Cdouble}}
    λ::Vector{Vector{Cdouble}}
    # t::Vector{Vector{Float64}}
end

function QuadraticProgramResult(stages::MultistageProblem)
    x = [zeros(stages.nx[i]) for i=1:stages.N+1]
    u = [zeros(stages.nu[i]) for i=1:stages.N+1]
    π = [zeros(stages.nx[i]) for i=2:stages.N+1]
    λ = [zeros(2*(stages.nb[i] + stages.nc[i])) for i=1:stages.N+1]
    return QuadraticProgramResult(x, u, π, λ)
end
