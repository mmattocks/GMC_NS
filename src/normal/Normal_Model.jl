mutable struct Normal_Record <: GMC_NS_Model_Record
    path::String
    log_Li::Float64
end

function construct_normal_record(path,log_Li)
    return Normal_Record(path, log_Li)
end

mutable struct Normal_Model <: GMC_NS_Model
    trajectory::Integer #id integer assigned at construction
    i::Integer #id of parent model, if any

    θ::Vector{Float64} #model's parameter Vector
    v::Vector{Float64} #model's velocity in parameter space
    log_Li::Float64 #model's likelihood, calculated by constructor

    function Normal_Model(trajectory::Integer, i::Integer, θ::Vector{Float64}, v::Vector{Float64}, obs::Vector{Float64}; v_init=false)
        θ=bound_σ!(θ)
        mod_normal=0
        try
            mod_normal=Normal(θ...)
        catch
            println(θ)
        end
        log_lh=lps(logpdf(mod_normal, obs))
        v_init && (v=rand(MvNormal(length(θ),1.)))

        new(trajectory, i, θ, v, log_lh)
    end
end

function construct_normal_model(args...; kwargs...)
    return Normal_Model(args...; kwargs...)
end

                function bound_σ!(θvec)
                    if θvec[2] <= 0.
                        θvec[2]=nextfloat(0.); return θvec
                    else
                        return θvec
                    end                        
                end


function Base.show(io::IO, m::Normal_Model; xsteps=100)
    μ,σ=m.θ
    n=Normal(μ,σ)
    X=[μ-4σ:((μ+4σ)-(μ-4σ))/xsteps:μ+4σ...]
    y=pdf.(n,X)

    show(io, lineplot(X,y,xlabel="X",ylabel="p",title="Normal Model $(m.id)"))
    println()
    println("Origin: $(m.origin)")
    println("v: $(m.v)")
end