mutable struct Eggbox_Ensemble <: GMC_NS_Ensemble
    path::String

    model_initλ::Function
    models::Vector{Eggbox_Record}

    contour::Float64
    log_Li::Vector{Float64}
    log_Xi::Vector{Float64}
    log_wi::Vector{Float64}
    log_Liwi::Vector{Float64}
    log_Zi::Vector{Float64}
    Hi::Vector{Float64}

    obs::Nothing
    priors::Vector{<:Distribution}
    constants::Vector{<:Real}
    box::Matrix{Float64}

    sample_posterior::Bool
    posterior_samples::Vector{Eggbox_Record}

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

Eggbox_Ensemble(path::String, no_models::Integer, prior, box, GS_settings...; sample_posterior::Bool=true) =
Eggbox_Ensemble(
    path,
    construct_eggbox_model,
    assemble_EMs(path, no_models, prior, box, Vector{Real}())...,
    [-Inf], #L0 = 0
	[0.], #X0 = 1
	[-Inf], #w0 = 0
	[-Inf], #Liwi0 = 0
	[-Inf], #Z0 = 0
	[0.], #H0 = 0,
    nothing,
    length(prior)==1 ? ([GMC_NS.marginals(prior)...]) : (prior),
    Vector{Real}(),
    length(prior)==1 ? (to_unit_ball.(box,[GMC_NS.marginals(prior)...])) : (to_unit_ball.(box,prior)),
    sample_posterior,
    Vector{Eggbox_Record}(),
    GS_settings...,
	no_models+1)

function Base.show(io::IO, e::Eggbox_Ensemble; progress=true)
    x1o=[0.0, 0.4, 0.8, 0.2, 0.6, 1.0, 0.0, 0.4, 0.8, 0.2, 0.6, 1.0, 0.0,0.4,0.8,0.2,0.6,1.0]
    x2o=[0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8,1.0, 1.0, 1.0]

    x1m=[m.pos[1] for m in e.models]
    x2m=[m.pos[2] for m in e.models]
    x1m.+=1.; x2m.+=1.
    x1m./=2.; x2m./=2

    plt=scatterplot(x1o,x2o,title="Eggbox Ensemble: Contour $(round(e.contour,digits=2))",color=:green,name="Optima")
    scatterplot!(plt, x1m, x2m, color=:magenta, name="Model Particles")
    show(io, plt)
    println()

    (progress && return 19)
end

function assemble_EMs(path::String, no_trajectories::Integer, prior, box, constants)
	ensemble_records = Vector{Eggbox_Record}()
    !isdir(path) && mkpath(path)
    length(prior)==1 ? (marginals=[GMC_NS.marginals(prior)...]) : (marginals=prior)
    box=to_unit_ball.(box,marginals)

    @showprogress 1 "Assembling Normal Model ensemble..." for trajectory_no in 1:no_trajectories
		model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            proposal=[rand.(prior)...]
            pos=to_unit_ball.(proposal,marginals)
            box_bound!(pos,box)
            θvec=to_prior.(pos,marginals)

            model = construct_eggbox_model(trajectory_no, 1, θvec, pos, [0.], nothing, constants...; v_init=true)
            
			serialize(model_path, model) #save the model to the ensemble directory
			push!(ensemble_records, Eggbox_Record(trajectory_no,1,pos,model_path,model.log_Li))
		else #interrupted assembly pick up from where we left off
			model = deserialize(model_path)
			push!(ensemble_records, Eggbox_Record(trajectory_no,1,model.pos,model_path,model.log_Li))
		end
	end

	return ensemble_records, minimum([record.log_Li for record in ensemble_records])
end