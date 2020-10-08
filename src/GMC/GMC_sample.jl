"""
    GMC_sample(e, τ)

Obtain a sample from the ensemble prior, constrained to the box defined by the current ensemble likelihood contour, by Galilean Monte Carlo, using ensemble particle velocities scaled by the time step 'τ', with 'r=1' (the least likely particle is discarded in favour of one inside its contour). Continue to generate new particles by GMC until one inside the contour is found or e.GMC_exhaust_σ iterates have failed to produce one. In the latter case, returns 'nothing'. In either case returns a BitVector history of calls to galilean_trajectory_sample! in which a 0 indicates a sampling failure and reversed particle.
"""
function GMC_sample(e::GMC_NS_Ensemble, τ::Float64)
    report=falses(0)
	for particle in 1:(Int64(floor(length(e.models)*e.GMC_exhaust_σ)))
		m_record = rand(e.models)
        m = deserialize(m_record.path)

        candidate=galilean_trajectory_sample!(m,e,τ)
        
        if candidate.id==m.id #particle reversed rather than new sample
            push!(report, false)
            serialize(m.path,candidate) #handle serialization here for a single process
        else
            push!(report, true)
            return candidate, report
        end
    end

    return nothing, report
end


function permute_IPM(e::IPM_Ensemble, job_chan::RemoteChannel, models_chan::RemoteChannel, comms_chan::RemoteChannel) #ensemble.models is partially updated on the worker to populate arguments for permute funcs
	persist=true
    id=myid()
    model_ctr=1
    put!(comms_chan,id)
    while persist
        wait(job_chan)
        start=time()
        e.models, e.contour, instruction = fetch(job_chan)
        instruction == "stop" && (persist=false; break)

        call_report=Vector{Tuple{Int64,Float64,Float64}}()
        for model=1:instruction.model_limit
			found::Bool=false
            m_record = rand(e.models)
            
            remotecall_fetch(isfile,1,m_record.path) ? m = (remotecall_fetch(deserialize,1,m_record.path)) : (break)

            model_blacklist=copy(m.permute_blacklist)
            filteridxs, filtered_funcs, filtered_weights, filtered_args = filter_permutes(instruction, model_blacklist)

            for call in 1:instruction.func_limit
                start=time()
                funcidx=rand(filtered_weights)
                permute_func=filtered_funcs[funcidx]
                pos_args,kw_args=get_permfunc_args(permute_func,e,m,filtered_args[funcidx])
                new_m=permute_func(pos_args...;kw_args...)
                push!(call_report,(filteridxs[funcidx],time()-start,new_m.log_Li - e.contour))

                if length(new_m.permute_blacklist) > 0 && new_m.log_Li == -Inf #in this case the blacklisted function will not work on this model at all; it should be removed from the functions to use
                    vcat(model_blacklist,new_m.permute_blacklist)
                    filteridxs, filtered_funcs, filtered_weights, filtered_args = filter_permutes(instruction, model_blacklist)
                end

                dupecheck(new_m,m) && new_m.log_Li > e.contour && ((put!(models_chan, (new_m ,id, call_report))); found=true; model_ctr=1; break)
                
			end
            found==true && break;
            model_ctr+=1
            wait(job_chan)
            fetch(job_chan)!=e.models && (break) #if the ensemble has changed during the search, update it
			model==instruction.model_limit && (put!(models_chan, ("quit", id, call_report));persist=false)#worker to put "quit" on channel if it fails to find a model more likely than contour
		end
	end
end
                function filter_permutes(instruction, model_blacklist)
                    filteridxs=findall(p->!in(p,model_blacklist),instruction.funcs)
                    filtered_funcs=instruction.funcs[filteridxs]
                    filtered_weights=filter_weights(instruction.weights, filteridxs)
                    filtered_args=instruction.args[filteridxs]
                    return filteridxs, filtered_funcs, filtered_weights, filtered_args
                end

                function filter_weights(weights, idxs)
                    deplete=0.
                    for (i, weight) in enumerate(weights)
                        !in(i, idxs) && (deplete+=weight)
                    end
                    weights[idxs].+=deplete/length(idxs)
                    new_weights=weights[idxs]./sum(weights[idxs])
                    return Categorical(new_weights)
                end

                function get_permfunc_args(func::Function,e::IPM_Ensemble, m::ICA_PWM_Model, args::Vector{Tuple{Symbol,Any}})
                    pos_args=[]
                    argparts=Base.arg_decl_parts(methods(func).ms[1])
                    argnames=[Symbol(argparts[2][n][1]) for n in 2:length(argparts[2])]
                    for argname in argnames #assemble basic positional arguments from ensemble and model fields
                        if argname == Symbol('m')
                            push!(pos_args,m)
                        elseif argname in fieldnames(IPM_Ensemble)
                            push!(pos_args,getfield(e,argname))
                        elseif argname in fieldnames(ICA_PWM_Model)
                            push!(pos_args,getfield(m,argname))
                        else
                            throw(ArgumentError("Positional argument $argname of $func not available in the ensemble or model!"))
                        end
                    end
            
                    kw_args = NamedTuple()
                    if length(args) > 0 #if there are any keyword arguments to pass
                        sym=Vector{Symbol}()
                        val=Vector{Any}()
                        for arg in args
                            push!(sym, arg[1]); push!(val, arg[2])
                        end
                        kw_args = (;zip(sym,val)...)
                    end

                    return pos_args, kw_args
                end

                function dupecheck(new_model, model)
                    (new_model.sources==model.sources && new_model.mix_matrix==model.mix_matrix) ? (return false) : (return true)
                end