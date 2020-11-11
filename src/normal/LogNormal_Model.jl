mutable struct LogNormal_Model <: GMC_NS_Model
    trajectory::Integer #id integer assigned at construction
    i::Integer #id of parent model, if any

    θ::Vector{Float64} #model's parameter Vector
    log_Li::Float64 #model's likelihood, calculated by constructor

    pos::Vector{Float64}
    v::Vector{Float64} #model's velocity in parameter space

    function LogNormal_Model(trajectory::Integer, i::Integer, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Float64}; v_init=false)
        μ,λ=θ
        mod_lnormal=LogNormal(μ,sqrt(1/λ))
        log_lh=lps(logpdf(mod_lnormal, obs))
        v_init && (v=rand(MvNormal(length(θ),1.)))

        new(trajectory, i, θ, log_lh, pos, v)
    end
end

function construct_lognormal_model(args...; kwargs...)
    return LogNormal_Model(args...; kwargs...)
end

function Base.show(io::IO, m::LogNormal_Model; xsteps=100)
    μ,λ=m.θ
    σ=sqrt(1/λ)
    n=LogNormal(μ,σ)
    X=[quantile(n,.025):(quantile(n,.975)-quantile(n,.025))/xsteps:quantile(n,.975)...]
    y=pdf.(n,X)

    show(io, lineplot(X,y,xlabel="X",ylabel="p",title=title="LogNormal Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)"))
    println()
    println("θ: $(m.θ)")
    println("v: $(m.v)")
end