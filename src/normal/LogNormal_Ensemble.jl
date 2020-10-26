mutable struct LogNormal_Ensemble <: GMC_NS_Ensemble
    path::String

    model_initλ::Function
    models::Vector{Normal_Record}

    contour::Float64
    log_Li::Vector{Float64}
    log_Xi::Vector{Float64}
    log_wi::Vector{Float64}
    log_Liwi::Vector{Float64}
    log_Zi::Vector{Float64}
    Hi::Vector{Float64}

    obs::Vector{Float64}
    priors::Vector{<:Distribution}
    constants::Vector{<:Real}

    sample_posterior::Bool
    posterior_samples::Vector{Normal_Record}

    GMC_tune_ι::Int64
    GMC_tune_ρ::Float64
    GMC_tune_μ::Int64
    GMC_tune_α::Float64
    GMC_tune_β::Float64

    GMC_timestep_η::Float64
    GMC_reflect_η::Float64
    GMC_exhaust_σ::Float64

    model_counter::Int64
end

LogNormal_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractFloat}, priors::AbstractVector{<:Distribution}, GMC_settings...; sample_posterior::Bool=true) =
LogNormal_Ensemble(
    path,
    construct_lognormal_model,
    assemble(construct_lognormal_model, construct_normal_record, path, no_models, obs, priors,Vector{Real}())...,
    [-Inf], #L0 = 0
	[0.], #ie exp(0) = all of the prior is covered
	[-Inf], #w0 = 0
	[-Inf], #Liwi0 = 0
	[-1e300], #Z0 = 0
	[0.], #H0 = 0,
    obs,
    priors,
    Vector{Real}(),
    sample_posterior,
    Vector{Normal_Record}(),
    GMC_settings...,
	no_models+1)

function Base.show(io::IO, m::LogNormal_Model, e::LogNormal_Ensemble; xsteps=100, progress=true)
    μ,σ=m.θ
    n=LogNormal(μ,σ)
    X=[quantile(n,.025):(quantile(n,.975)-quantile(n,.025))/xsteps:quantile(n,.975)...]
    y=pdf.(n,X)

    scattermin=maximum(y)*.25
    scattermax=maximum(y)*.75
    scatterrange=[scattermin:(scattermax-scattermin)/length(e.obs):scattermax...]
    scattery=[rand(scatterrange) for i in 1:length(e.obs)]

    xlims=(min(minimum(e.obs),minimum(X)),max(maximum(e.obs),maximum(X)))

    plt=lineplot(X,y,xlabel="X",ylabel="p",title="LogNormal Model $(m.trajectory).$(m.i), log_Li $(m.log_Li)", name="Model", xlim=xlims)
    scatterplot!(plt, e.obs, scattery, name="Obs")

    show(io, plt)
    println() 
    println("θ: $(m.θ)")
    println("v: $(m.v)")

    progress && return nrows(plt.graphics)+8
end