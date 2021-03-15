function galilean_search!(m,e,t,llidx,samples,s_limit,cache)
    new_m,cache=galilean_trajectory_sample!(m,e,t,cache)
    m.i==1 && ((m,new_m,cache)=start_trim(m,new_m,cache))


    while new_m.i==m.i && samples<=s_limit && t.τ > e.GMC_τ_death#failure to find any new step, retune tau and try again until sample limit reached or trajectory's tuner reaches taudeath
        new_m,cache=galilean_trajectory_sample!(m,e,t,cache)
        samples+=1
        m.i==1 && ((m,new_m,cache)=start_trim(m,new_m,cache))
    end

    if new_m.i!=m.i
        e.sample_posterior && push!(e.posterior_samples, e.models[llidx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
        rectype=typeof(e.models[1])
        deleteat!(e.models, llidx) #remove least likely model record

        new_model_record = rectype(new_m.trajectory,new_m.i, new_m.pos, string(e.path,'/',new_m.trajectory,'.',new_m.i), new_m.log_Li)
        push!(e.models, new_model_record) #insert new model record
        serialize(new_model_record.path, new_m)
        
        return true, cache
    else
        return false, cache
    end
end

function start_trim(m,new_m,cache)
    if new_m.i!=m.i && isapprox(new_m.pos,m.pos)
        params=[getfield(m,pn) for pn in propertynames(m)]
        params[6]=new_m.v
        m=typeof(m)(params...)
        params=[getfield(cache,pn) for pn in propertynames(m)]
        params[2]=2;
        new_m=typeof(m)(params...)
        cache=nothing   
    end
    return m,new_m,cache
end

function diffusion_search!(e, tuners, llidx, sample_limit, τ_limit)
    sample=0; model_selected=false; new_m=0; src_traj=0
    while !(sample > sample_limit) && !model_selected
        sample +=1
        mrec=rand(e.models)
        model=deserialize(mrec.path)
        randdir=rand(MvNormal(length(model.θ),1.))
        τ=tuners[model.trajectory].τ
        new_pos=model.pos+randdir*τ
        box_bound!(new_pos,e.box)
        new_m=e.model_initλ(e.t_counter, 1, to_prior.(new_pos,e.priors), new_pos, randdir, e.obs, e.constants...)

        while τ*.7>τ_limit && new_m.log_Li<e.contour
            τ*=.7
            randdir=rand(MvNormal(length(model.θ),1.))
            new_pos=model.pos+randdir*τ
            box_bound!(new_pos,e.box)
            new_m=e.model_initλ(e.t_counter, 1, to_prior.(new_pos,e.priors), new_pos, randdir, e.obs, e.constants...)
        end

        new_m.log_Li >= e.contour && (model_selected=true; src_traj=model.trajectory)
    end

    if new_m.log_Li >= e.contour
        e.sample_posterior && push!(e.posterior_samples, e.models[llidx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
        rectype=typeof(e.models[1])
        deleteat!(e.models, llidx) #remove least likely model record
        
        new_model_record = rectype(new_m.trajectory,new_m.i, new_m.pos, string(e.path,'/',new_m.trajectory,'.',new_m.i), new_m.log_Li)
        push!(e.models, new_model_record) #insert new model record
        serialize(new_model_record.path, new_m)
        tuners[e.t_counter]=τ_PID(e)
        tuners[e.t_counter].τ=tuners[src_traj].τ
        e.t_counter +=1
        return true
    else
        return false
    end
end

function insert_reflection!(e, cache, llidx)
    e.sample_posterior && push!(e.posterior_samples, e.models[llidx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
    rectype=typeof(e.models[1])
    deleteat!(e.models, llidx) #remove least likely model record

    new_model_record = rectype(cache.trajectory,cache.i, cache.pos, string(e.path,'/',cache.trajectory,'.',cache.i), cache.log_Li)
    push!(e.models, new_model_record) #insert new model record
    serialize(new_model_record.path, cache)
end

function ellipsoid_resample!(e,tdict,llidx,samples,s_limit)
    posmat=zeros(length(e.models[1].pos),length(e.models))
    for (n,rec) in enumerate(e.models)
        posmat[:,n]=rec.pos
    end
    e_ellip=bounding_ellipsoid(posmat)

    new_pos=sample_ellipsoid(e_ellip)
    box_bound!(new_pos,e.box)
    new_m=e.model_initλ(e.t_counter, 1, to_prior.(new_pos,e.priors), new_pos, [0.], e.obs, e.box, e.constants..., v_init=true)
    samples+=1

    while new_m.log_Li < e.contour && samples < s_limit
        new_pos=sample_ellipsoid(e_ellip)
        box_bound!(new_pos,e.box)
        new_m=e.model_initλ(e.t_counter, 1, to_prior.(new_pos,e.priors), new_pos, [0.], e.obs, e.box, e.constants..., v_init=true)
        samples+=1
    end

    if new_m.log_Li > e.contour
        e.sample_posterior && push!(e.posterior_samples, e.models[llidx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
        rectype=typeof(e.models[1])
        deleteat!(e.models, llidx) #remove least likely model record
        
        new_model_record = rectype(new_m.trajectory,new_m.i, new_m.pos, string(e.path,'/',new_m.trajectory,'.',new_m.i), new_m.log_Li)
        push!(e.models, new_model_record) #insert new model record
        serialize(new_model_record.path, new_m)
        tdict[e.t_counter]=τ_PID(e)
        e.t_counter +=1
        return true
    else
        return false
    end
end