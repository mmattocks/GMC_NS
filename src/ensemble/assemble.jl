function assemble(mλ::Function, rλ::Function, path::String, no_trajectories::Integer, obs, priors, constants)
    rectest=rλ("test",0.)
	ensemble_records = Vector{typeof(rλ("test",0.))}()
	!isdir(path) && mkpath(path)
    @showprogress 1 "Assembling $mλ ensemble..." for trajectory_no in 1:no_trajectories
		model_path = string(path,'/',trajectory_no,'.',1)
        if !isfile(model_path)
            θvec=sample_θvec(priors)
            model = mλ(trajectory_no, 1, θvec, [0.], obs, constants...; v_init=true)
            while model.log_Li == -Inf
                θvec=sample_θvec(priors)
                model = mλ(trajectory_no, 1, θvec, [0.], obs, constants...; v_init=true)
            end
			serialize(model_path, model) #save the model to the ensemble directory
			push!(ensemble_records, rλ(model_path,model.log_Li))
		else #interrupted assembly pick up from where we left off
			model = deserialize(model_path)
			push!(ensemble_records, rλ(model_path,model.log_Li))
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