type QuadraticProgram
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

    QuadraticProgram(N::Int) = new(N)
end

function QuadraticProgram(N::Int, nx::Int, nu::Int, nb::Int, nc::Int)
    nxvect = fill(nx, N+1)
    nuvect = fill(nu, N+1)
    nbvect = fill(nb, N+1)
    ncvect = fill(nc, N+1)
    nuvect[N+1] = 0
    if nb == 0; nbvect[1] = nx; end
    qp = QuadraticProgram(N, nxvect, nuvect, nbvect, ncvect)
    if nb == 0; qp.idx[1] .= nu:nx+nu-1; end
    return qp
end

function QuadraticProgram(N::Int, nx::Vector{Int}, nu::Vector{Int},
                          nb::Vector{Int}, nc::Vector{Int})
    length(nx) == N+1 || error()
    length(nu) == N+1 || error()
    length(nb) == N+1 || error()
    length(nc) == N+1 || error()

    qp = QuadraticProgram(N)

    qp.nx = nx
    qp.nu = nu
    qp.nb = nb
    qp.nc = nc

    qp.A = [zeros(nx[i], nx[i]) for i=1:N]
    qp.B = [zeros(nx[i], nu[i]) for i=1:N]
    qp.b = [zeros(nx[i]) for i=1:N]

    qp.Q = [zeros(nx[i], nx[i]) for i=1:N+1]
    qp.S = [zeros(nu[i], nx[i]) for i=1:N]
    qp.R = [zeros(nu[i], nu[i]) for i=1:N]
    qp.q = [zeros(nx[i]) for i=1:N+1]
    qp.r = [zeros(nu[i]) for i=1:N]

    qp.lb = [zeros(nb[i]) for i=1:N+1]
    qp.ub = [zeros(nb[i]) for i=1:N+1]
    qp.idx = [zeros(Int, nb[i]) for i=1:N+1]

    qp.lc = [zeros(nc[i]) for i=1:N+1]
    qp.uc = [zeros(nc[i]) for i=1:N+1]
    qp.C = [zeros(nc[i], nx[i]) for i=1:N+1]
    qp.D = [zeros(nc[i], nu[i]) for i=1:N]

    return qp
end

type QuadraticProgramResult
    x::Vector{Vector{Cdouble}}
    u::Vector{Vector{Cdouble}}
    π::Vector{Vector{Cdouble}}
    λ::Vector{Vector{Cdouble}}
    # t::Vector{Vector{Float64}}
end

function QuadraticProgramResult(qp::QuadraticProgram)
    x = [zeros(qp.nx[i]) for i=1:qp.N+1]
    u = [zeros(qp.nu[i]) for i=1:qp.N+1]
    π = [zeros(qp.nx[i]) for i=2:qp.N+1]
    λ = [zeros(2*(qp.nb[i] + qp.nc[i])) for i=1:qp.N+1]
    return QuadraticProgramResult(x, u, π, λ)
end
