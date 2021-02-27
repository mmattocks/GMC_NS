mutable struct Normal_Record <: GMC_NS_Model_Record
    trajectory::Int64
    i::Int64
    pos::Vector{Float64}
    path::String
    log_Li::Float64
end

struct Normal_Model <: GMC_NS_Model
    trajectory::Integer #id integer assigned at construction
    i::Integer #id of parent model, if any

    θ::Vector{Float64} #model's parameter Vector
    log_Li::Float64 #model's likelihood, calculated by constructor

    pos::Vector{Float64}
    v::Vector{Float64} #model's velocity in parameter space
end

function construct_normal_model(trajectory::Integer, i::Integer, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Float64}; v_init=false)
    μ,λ=θ
    mod_normal=Normal(μ,sqrt(1/λ))
    log_lh=lps(logpdf.(mod_normal, obs))
    v_init && (v=rand(MvNormal(length(θ),1.)))

    Normal_Model(trajectory, i, θ, log_lh, pos, v)
end

function Base.show(io::IO, m::Normal_Model; xsteps=100)
    μ,λ=m.θ
    σ=sqrt(1/λ)
    n=Normal(μ,σ)
    X=[μ-4σ:((μ+4σ)-(μ-4σ))/xsteps:μ+4σ...]
    y=pdf.(n,X)

    show(io, lineplot(X,y,xlabel="X",ylabel="p",title="Normal Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)"))
    println()
    println("θ: $(m.θ)")
    println("v: $(m.v)")
end