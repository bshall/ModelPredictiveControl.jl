type HPMPCSolver <: AbstractSolver
    tol::Float64
    maxit::Int
    μ₀::Float64
    warmstart::Bool
    # ux₀::Vector{Vector{Float64}}
    # π₀::Vector{Vector{Float64}}
    # λ₀::Vector{Vector{Float64}}
    # t₀::Vector{Vector{Float64}}
    res::Vector{Float64}
    stats::Matrix{Float64}
    workspace::Vector{Float64}

    function HPMPCSolver(workspace::AbstractVector, tol::Real=1e-8, maxit::Int=10,
                         μ₀::Real=0.0, warmstart::Bool=false)
        res = zeros(5)
        stats = zeros(maxit, 5)
        return new(tol, maxit, μ₀, warmstart, res, stats, workspace)
    end
end

function getworkspacesize(::Type{HPMPCSolver}, qp::QuadraticProgram)
    size = ccall((:hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes, libhpmpc),
                 Cint, (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                 Ptr{Ptr{Cint}}, Ptr{Cint}, Cint),
                 qp.N, qp.nx, qp.nu, qp.nb, qp.idx, qp.nc, qp.N)
end
