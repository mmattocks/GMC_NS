mutable struct Normal_Ensemble <: GMC_NS_Ensemble
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

Normal_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractFloat}, priors::AbstractVector{<:Distribution}, GMC_settings...; sample_posterior::Bool=true) =
Normal_Ensemble(
    path,
    construct_normal,
    assemble_NMs(path, no_models, obs, priors)...,
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

function assemble_NMs(path::String, no_models::Integer, obs, priors)
	ensemble_records = Vector{Normal_Record}()
	!isdir(path) && mkpath(path)

    for model_no in 1:no_models
		model_path = string(path,'/',model_no)
        if !isfile(model_path)
            θvec=sample_θvec(priors)
			model = Normal_Model(model_no, 0, θvec, [0.], obs; v_init=true)
			serialize(model_path, model) #save the model to the ensemble directory
			push!(ensemble_records, Normal_Record(model_path,model.log_Li))
		else #interrupted assembly pick up from where we left off
			model = deserialize(model_path)
			push!(ensemble_records, Normal_Record(model_path,model.log_Li))
		end
	end

	return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end

                function sample_θvec(priors)
                    θvec=zeros(0)
                    for prior in priors
                        sample=rand(prior)
                        for val in sample
                            push!(θvec,val)
                        end
                    end
                    return θvec
                end

function Base.show(io::IO, m::Normal_Model, e::Normal_Ensemble; xsteps=100, progress=true)
    μ,σ=m.θ
    n=Normal(μ,σ)
    X=[μ-4σ:((μ+4σ)-(μ-4σ))/xsteps:μ+4σ...]
    y=pdf.(n,X)

    scattermin=y[Int(floor(xsteps/2)-floor(xsteps/3))]
    scattermax=y[Int(floor(xsteps/2)-floor(xsteps/6))]
    scatterrange=[scattermin:(scattermax-scattermin)/length(e.obs):scattermax...]
    scattery=[rand(scatterrange) for i in 1:length(e.obs)]

    xlims=(min(minimum(e.obs),minimum(X)),max(maximum(e.obs),maximum(X)))

    plt=lineplot(X,y,xlabel="X",ylabel="p",title="Normal Model $(m.id)", name="Model", xlim=xlims)
    scatterplot!(plt, e.obs, scattery, name="Obs")

    show(io, plt)
    println()
    println("Origin: $(m.origin)")
    println("v: $(m.v)")

    progress && return nrows(plt.graphics)+8
end