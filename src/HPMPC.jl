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

function getworkspacesize(::Type{HPMPCSolver}, stages::MultistageProblem)
    size = ccall((:hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes, libhpmpc),
                 Cint, (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                 Ptr{Ptr{Cint}}, Ptr{Cint}, Cint),
                 stages.N, stages.nx, stages.nu, stages.nb,
                 stages.idx, stages.nc, stages.N)
end

function solve!(out::QuadraticProgramResult, stages::MultistageProblem,
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
                   stages.N, stages.nx, stages.nu, stages.nb,
                   stages.idx, stages.nc, stages.N, solver.warmstart,
                   stages.A, stages.B,
                   stages.b, stages.Q,
                   stages.S, stages.R,
                   stages.q, stages.r,
                   stages.lb, stages.ub,
                   stages.C, stages.D,
                   stages.lc, stages.uc,
                   out.x, out.u,
                   out.π, out.λ,
                   solver.res, solver.workspace, solver.stats)
end
