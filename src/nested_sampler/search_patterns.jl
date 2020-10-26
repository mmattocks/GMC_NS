function ellipsoid_resample!(e,tdict,llidx,samples,s_limit)
    posmat=zeros(length(e.models[1].pos),length(e.models))
    for (n,rec) in enumerate(e.models)
        posmat[:,n]=rec.pos
    end
    e_ellip=bounding_ellipsoid(posmat)

    @info "Ellipsoid resample $(samples)"

    new_pos=sample_ellipsoid(e_ellip)
    box_bound!(new_pos,e.box)
    new_m=e.model_initλ(e.t_counter, 1, to_prior.(new_pos,e.priors), new_pos, [0.], e.obs, e.box, e.constants..., v_init=true)
    samples+=1

    while new_m.log_Li < e.contour && samples < s_limit
        move_cursor_up_while_clearing_lines(stdout,1)
        @info "Ellipsoid resample $(samples)"
        new_pos=sample_ellipsoid(e_ellip)
        box_bound!(new_pos,e.box)
        new_m=e.model_initλ(e.t_counter, 1, to_prior.(new_pos,e.priors), new_pos, [0.], e.obs, e.box, e.constants..., v_init=true)
        samples+=1
    end

    if new_m.log_Li > e.contour
        move_cursor_up_while_clearing_lines(stdout,1)
        @info "Found"
        e.sample_posterior && push!(e.posterior_samples, e.models[llidx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
        deleteat!(e.models, llidx) #remove least likely model record

        new_model_record = typeof(e.models[1])(new_m.trajectory,new_m.i, new_m.pos, string(e.path,'/',new_m.trajectory,'.',new_m.i), new_m.log_Li)
        push!(e.models, new_model_record) #insert new model record
        serialize(new_model_record.path, new_m)
        tdict[e.t_counter]=τ_PID(e)
        e.t_counter +=1
        move_cursor_up_while_clearing_lines(stdout,1)
        return true
    else
        move_cursor_up_while_clearing_lines(stdout,1)
        @warn "Ellipsoid resample failure!"
        move_cursor_up_while_clearing_lines(stdout,1)
        return false
    end
end


function galilean_search!(m,e,t,llidx,samples,s_limit)
    candidate=galilean_trajectory_sample!(m,e,t)

    @info "Galilean search $(samples), trajectory $(m.trajectory), i $(m.i), τ $(t.τ) † $(e.GMC_τ_death)"

    while candidate.i==m.i && samples<=s_limit && t.τ > e.GMC_τ_death#failure to find any new step, retune tau and try again until sample limit reached
        move_cursor_up_while_clearing_lines(stdout,1)
        @info "Galilean search $(samples), trajectory $(m.trajectory), i $(m.i), τ $(t.τ) † $(e.GMC_τ_death)"

        candidate=galilean_trajectory_sample!(m,e,t)
        samples+=1
        samples%4==0 && (m.v=rand(MvNormal(length(m.θ),1.)))
    end

    if candidate.i!=m.i
        move_cursor_up_while_clearing_lines(stdout,1)
        @info "Found"

        e.sample_posterior && push!(e.posterior_samples, e.models[llidx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
        deleteat!(e.models, llidx) #remove least likely model record

        new_model_record = typeof(e.models[1])(candidate.trajectory,candidate.i, candidate.pos, string(e.path,'/',candidate.trajectory,'.',candidate.i), candidate.log_Li)
        push!(e.models, new_model_record) #insert new model record
        serialize(new_model_record.path, candidate)
        move_cursor_up_while_clearing_lines(stdout,1)

        return true
    else
        move_cursor_up_while_clearing_lines(stdout,1)
        @info "†"
        move_cursor_up_while_clearing_lines(stdout,1)
        return false
    end
end