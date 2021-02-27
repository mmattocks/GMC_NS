#NOTE: LogNormal_Models reuse Normal_Records

struct LogNormal_Model <: GMC_NS_Model
    trajectory::Integer #id integer assigned at construction
    i::Integer #id of parent model, if any

    θ::Vector{Float64} #model's parameter Vector
    log_Li::Float64 #model's likelihood, calculated by constructor

    pos::Vector{Float64}
    v::Vector{Float64} #model's velocity in parameter space
end

function construct_lognormal_model(trajectory::Integer, i::Integer, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Float64}; v_init=false)
    μ,λ=θ
    lμ=log(μ^2/sqrt(μ^2 + inv(λ)))
    lσ=sqrt(log(1+(inv(λ)/μ^2)))
    mod_lnormal=LogNormal(lμ,lσ)
    log_lh=lps(logpdf.(mod_lnormal, obs))
    v_init && (v=rand(MvNormal(length(θ),1.)))

    LogNormal_Model(trajectory, i, θ, log_lh, pos, v)
end

function Base.show(io::IO, m::LogNormal_Model; xsteps=100)
    μ,λ=m.θ
    lμ=log(μ^2/sqrt(μ^2 + inv(λ)))
    lσ=sqrt(log(1+(inv(λ)/μ^2)))
    n=LogNormal(lμ,lσ)
    X=[quantile(n,.025):(quantile(n,.975)-quantile(n,.025))/xsteps:quantile(n,.975)...]
    y=pdf.(n,X)

    show(io, lineplot(X,y,xlabel="X",ylabel="p",title=title="LogNormal Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)"))
    println()
    println("θ: $(m.θ)")
    println("v: $(m.v)")
end