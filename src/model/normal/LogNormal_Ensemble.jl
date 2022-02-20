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
    box::Matrix{Float64}

    sample_posterior::Bool
    posterior_samples::Vector{Normal_Record}

    GMC_Nmin::Int64

    GMC_τ_death::Float64
    GMC_init_τ::Float64
    GMC_tune_μ::Int64
    GMC_tune_α::Float64
    GMC_tune_PID::NTuple{3,Float64}

    GMC_timestep_η::Float64
    GMC_reflect_η::Float64
    GMC_exhaust_σ::Float64
    GMC_chain_κ::Int64

    t_counter::Int64
end

LogNormal_Ensemble(path::String, no_models::Integer, obs::AbstractVector{<:AbstractFloat}, prior, box, GS_settings...; sample_posterior::Bool=true) =
LogNormal_Ensemble(
    path,
    construct_lognormal_model,
    assemble_lNMs(path, no_models, obs, prior, box, Vector{Real}())...,
    [-Inf], #L0 = 0
	[0.], #X0 = 1
	[-Inf], #w0 = 0
	[-Inf], #Liwi0 = 0
	[-Inf], #Z0 = 0
	[0.], #H0 = 0,
    obs,
    length(prior)==1 ? ([GMC_NS.marginals(prior)...]) : (prior),
    Vector{Real}(),
    length(prior)==1 ? (to_unit_ball.(box,[GMC_NS.marginals(prior)...])) : (to_unit_ball.(box,prior)),
    sample_posterior,
    Vector{Normal_Record}(),
    GS_settings...,
	no_models+1)

function Base.show(io::IO, m::LogNormal_Model, e::LogNormal_Ensemble; xsteps=100, progress=true)
    μ,λ=m.θ
    lμ=log(μ^2/sqrt(μ^2 + inv(λ)))
    lσ=sqrt(log(1+(inv(λ)/μ^2)))
    n=LogNormal(lμ,lσ)
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

function assemble_lNMs(path::String, no_trajectories::Integer, obs, prior, box, constants)
	ensemble_records = Vector{Normal_Record}()
    !isdir(path) && mkpath(path)
    length(prior)==1 ? (marginals=[GMC_NS.marginals(prior)...]) : (marginals=prior)
    box=to_unit_ball.(box,marginals)

    @showprogress 1 "Assembling log-Normal Model ensemble..." for trajectory_no in 1:no_trajectories
		model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=[rand.(prior)...]
            pos=to_unit_ball.(proposal,marginals)
            box_bound!(pos,box)
            θvec=to_prior.(pos,marginals)

            model = construct_lognormal_model(trajectory_no, 1, θvec, pos, [0.], obs, constants...; v_init=true)
            
			serialize(model_path, model) #save the model to the ensemble directory
			push!(ensemble_records, Normal_Record(trajectory_no,1,pos,model_path,model.log_Li))
		else #interrupted assembly pick up from where we left off
			model = deserialize(model_path)
			push!(ensemble_records, Normal_Record(trajectory_no,1,model.pos,model_path,model.log_Li))
		end
	end

	return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end