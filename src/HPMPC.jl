type HPMPCSolver <: AbstractSolver
    tol::Cdouble
    maxit::Cint
    μ₀::Cdouble
    warmstart::Bool
    # ux₀::Vector{Vector{Float64}}
    # π₀::Vector{Vector{Float64}}
    # λ₀::Vector{Vector{Float64}}
    # t₀::Vector{Vector{Float64}}
    res::Vector{Cdouble}
    stats::Matrix{Cdouble}
    workspace::Vector{Cdouble}

    function HPMPCSolver(workspace::AbstractVector, tol::Real=1e-8, maxit::Int=10,
                         μ₀::Real=0.0, warmstart::Bool=false)
        res = zeros(5)
        stats = zeros(5, maxit)
        return new(tol, maxit, μ₀, warmstart, res, stats, workspace)
    end
end

function getworkspacesize(::Type{HPMPCSolver}, qp::QuadraticProgram)
    size = ccall((:hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes, libhpmpc),
                 Cint, (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                 Ptr{Ptr{Cint}}, Ptr{Cint}, Cint),
                 qp.N, qp.nx, qp.nu, qp.nb, qp.idx, qp.nc, qp.N)
end

function solve!(out::QuadraticProgramResult, qp::QuadraticProgram,
                solver::HPMPCSolver)

    nit = Ref{Cint}(0)
    status = ccall((:fortran_order_d_ip_ocp_hard_tv, libhpmpc), Cint,
                   (Ptr{Cint}, Cint, Cdouble, Cdouble,
                    Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                    Ptr{Ptr{Cint}}, Ptr{Cint}, Cint, Cint,
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}},
                    Ptr{Cdouble}, Ptr{Void}, Ptr{Cdouble}),
                   nit, solver.maxit, solver.μ₀, solver.tol,
                   qp.N, qp.nx, qp.nu, qp.nb,
                   qp.idx, qp.nc, qp.N, solver.warmstart,
                   qp.A, qp.B,
                   qp.b, qp.Q,
                   qp.S, qp.R,
                   qp.q, qp.r,
                   qp.lb, qp.ub,
                   qp.C, qp.D,
                   qp.lc, qp.uc,
                   out.x, out.u,
                   out.π, out.λ,
                   solver.res, solver.workspace, solver.stats)
end
