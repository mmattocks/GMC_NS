mutable struct Eggbox_Record <: GMC_NS_Model_Record
    trajectory::Int64
    i::Int64
    pos::Vector{Float64}
    path::String
    log_Li::Float64
end

struct Eggbox_Model <: GMC_NS_Model
    trajectory::Integer #id integer assigned at construction
    i::Integer #id of parent model, if any

    θ::Vector{Float64} #model's parameter Vector
    log_Li::Float64 #model's likelihood, calculated by constructor

    pos::Vector{Float64}
    v::Vector{Float64} #model's velocity in parameter space
end

function construct_eggbox_model(trajectory::Integer, i::Integer, θ::Vector{Float64}, pos::Vector{Float64}, v::Vector{Float64}, obs::Vector{Float64}; v_init=false)
    x1,x2=θ
    log_lh=(2 + cos(5π * x1) * cos(5π * x2))^5
    v_init && (v=rand(MvNormal(length(θ),1.)))

    Eggbox_Model(trajectory, i, θ, log_lh, pos, v)
end

function Base.show(io::IO, m::Eggbox_Model; xsteps=100)
    println("Eggbox Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)")
    println("θ: $(m.θ)")
    println("v: $(m.v)")
end