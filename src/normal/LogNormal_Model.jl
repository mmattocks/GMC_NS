mutable struct LogNormal_Model <: GMC_NS_Model
    trajectory::Integer #id integer assigned at construction
    i::Integer #id of parent model, if any

    θ::Vector{Float64} #model's parameter Vector
    v::Vector{Float64} #model's velocity in parameter space
    log_Li::Float64 #model's likelihood, calculated by constructor

    function LogNormal_Model(trajectory::Integer, i::Integer, θ::Vector{Float64}, v::Vector{Float64}, obs::Vector{Float64}; v_init=false)
        bound_lnθ!(θ)
        mod_lnormal=0
        try
            mod_lnormal=LogNormal(θ...)
        catch
            println(θ)
        end
        log_lh=lps(logpdf(mod_lnormal, obs))
        v_init && (v=rand(MvNormal(length(θ),1.)))

        new(trajectory, i, θ, v, log_lh)
    end
end

function construct_lognormal_model(args...; kwargs...)
    return LogNormal_Model(args...; kwargs...)
end

                function bound_lnθ!(θ)
                    θ[findall(θi->θi<0., θ)].=nextfloat(0.)
                end


function Base.show(io::IO, m::LogNormal_Model; xsteps=100)
    μ,σ=m.θ
    n=LogNormal(μ,σ)
    X=[quantile(n,.025):(quantile(n,.975)-quantile(n,.025))/xsteps:quantile(n,.975)...]
    y=pdf.(n,X)

    show(io, lineplot(X,y,xlabel="X",ylabel="p",title=title="LogNormal Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)"))
    println()
        println("θ: $(m.θ)")
    println("v: $(m.v)")

end