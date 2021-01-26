function ensemble_history(e::GMC_NS_Ensemble, bins=25)
    !e.sample_posterior && throw(ArgumentError("This ensemble has no posterior samples to show a history for!"))
    livec=vcat([model.log_Li for model in e.models],[model.log_Li for model in e.posterior_samples])
    show(histogram(livec, nbins=bins))
end

function e_backup(e::GMC_NS_Ensemble, tuner::Dict{Int64,<:GMC_Tuner})
    serialize(string(e.path,'/',"ens"), e)
    serialize(string(e.path,'/',"tuner"), tuner)
end

function clean_ensemble_dir!(e::GMC_NS_Ensemble)
    if e.sample_posterior
        ens_paths=vcat(
            [basename(m.path) for m in e.models],
            [basename(m.path) for m in e.posterior_samples],
            "ens", "tuner"
        )
    else
        ens_paths=vcat(
            [basename(m.path) for m in e.models],
            "ens", "tuner"
        )
    end

    for file in readdir(e.path)
        !(file in ens_paths) && rm(e.path*'/'*file)
    end
end

function complete_evidence(e::GMC_NS_Ensemble)
    return final_logZ = logaddexp(e.log_Zi[end], (logsumexp([model.log_Li for model in e.models] .+  (e.log_Xi[length(e.log_Li)] - log(length(e.models))))))
end

function measure_evidence(e::GMC_NS_Ensemble)
    return ms=measurement(complete_evidence(e),sqrt(abs(e.Hi[end])/length(e.models)))
end


function reset_ensemble!(e::GMC_NS_Ensemble)
    Ni=Int64.(round.(collect(-1:-1:-length(e.log_Li)+1)./e.log_Xi[2:end]))
    N=Ni[1]

    new_models=Vector{typeof(e.models[1])}()

    for traj in 1:N
        if string(traj)*".1" in [basename(record.path) for record in e.models]
            push!(new_models,e.models[findfirst(isequal(string(traj)*".1"), [basename(record.path) for record in e.models])])
        else
            push!(new_models,e.posterior_samples[findfirst(isequal(string(traj)*".1"), [basename(record.path) for record in e.posterior_samples])])
        end
    end

    e.models=new_models
    e.posterior_samples=Vector{typeof(e.models[1])}()

    e.contour=minimum([record.log_Li for record in e.models])

    e.log_Li=[e.log_Li[1]]
    e.log_Xi=[e.log_Xi[1]]
    e.log_wi=[e.log_wi[1]]
    e.log_Liwi=[e.log_Liwi[1]]
    e.log_Zi=[e.log_Zi[1]]
    e.Hi=[e.Hi[1]]

    e.t_counter=length(e.models)+1

    clean_ensemble_dir(e, 0; ignore_warn=true)
    isfile(e.path*"/tuner") && rm(e.path*"/tuner")
    serialize(e.path*"/ens", e)

    return e
end

function move_ensemble!(e::GMC_NS_Ensemble,path::String)
    !isdir(path) && mkpath(path)
    rectype=typeof(e.models[1])

    for file in readdir(e.path)
        mv(e.path*'/'*file,path*'/'*file)
    end

    for (n,model) in enumerate(e.models)
        e.models[n]=rectype(model.trajectory, model.i, model.pos, path*'/'*basename(model.path), model.log_Li)
    end
    if e.sample_posterior
        for (n,model) in enumerate(e.posterior_samples)
            e.posterior_samples[n]=rectype(model.trajectory, model.i, model.pos, path*'/'*basename(model.path), model.log_Li)
        end
    end

    rm(e.path)
    e.path=path
    serialize(e.path*"/ens",e)
    return e
end

function copy_ensemble(e::GMC_NS_Ensemble,path::String)
    new_e=deepcopy(e)
    !isdir(path) && mkdir(path)
    rectype=typeof(e.models[1])
    for file in readdir(e.path)
        cp(e.path*'/'*file,path*'/'*file)
    end

    for (n,model) in enumerate(e.models)
        new_e.models[n]=rectype(model.trajectory, model.i, model.pos, path*'/'*basename(model.path), model.log_Li)
    end
    if e.sample_posterior
        for (n,model) in enumerate(e.posterior_samples)
            new_e.posterior_samples[n]=rectype(model.trajectory, model.i, model.pos, path*'/'*basename(model.path), model.log_Li)
        end
    end

    new_e.path=path
    serialize(new_e.path*"/ens",new_e)
    return new_e
end

function rewind_ensemble!(e::GMC_NS_Ensemble,rewind_idx)
    !e.sample_posterior && throw(ArgumentError("An ensemble not retaining posterior samples cannot be rewound!"))
    rewind_idx >= length(e.log_Li) && throw(ArgumentError("rewind_idx must be less than the current iterate!"))

    Ni=Int64.(round.(collect(-1:-1:-length(e.log_Li)+1)./e.log_Xi[2:end]))
    insert!(Ni,1,Ni[1])

    for it in length(e.log_Li)-1:-1:rewind_idx
        last_post=pop!(e.posterior_samples)

        if Ni[it]==Ni[it-1]
            models_idx=findfirst(isequal(last_post.trajectory), [rec.trajectory for rec in e.models])
            insert!(e.models, models_idx, last_post)
            deleteat!(e.models, models_idx+1)
        else
            push!(e.models, last_post)
        end

        !(length(e.models)==Ni[it-1]) && throw(DomainError("it: $it Ni[it]:$(Ni[it]) Ni[it-1]:$(Ni[it-1])"))
    end

    e.contour=e.log_Li[rewind_idx]
    e.log_Li=e.log_Li[1:rewind_idx]
    e.log_Xi=e.log_Xi[1:rewind_idx]
    e.log_wi=e.log_wi[1:rewind_idx]
    e.log_Liwi=e.log_Liwi[1:rewind_idx]
    e.log_Zi=e.log_Zi[1:rewind_idx]
    e.Hi=e.Hi[1:rewind_idx]

    traj_nos=vcat([rec.trajectory for rec in e.models],[rec.trajectory for rec in e.posterior_samples])

    e.t_counter=maximum(traj_nos)

    return e
end

function show_models(e::GMC_NS_Ensemble,idxs)
    liperm=sortperm([model.log_Li for model in e.models],rev=true)
    for idx in idxs
        m=deserialize(e.models[liperm[idx]].path)
        show(m)
    end
end

function show_models_e(e::GMC_NS_Ensemble,idxs)
    liperm=sortperm([model.log_Li for model in e.models],rev=true)
    for idx in idxs
        m=deserialize(e.models[liperm[idx]].path)
        show(stdout,m,e)
    end
end

function get_model(e::GMC_NS_Ensemble,no)
    return deserialize(e.path*'/'*string(no))
end

#recalculate the ensemble's history from the succession of log_Li and log_Xi values
function recalculate_ensemble!(e::GMC_NS_Ensemble, wi_mode="trapezoidal")
    Ni=Int64.(round.(collect(-1:-1:-length(e.log_Li)+1)./e.log_Xi[2:end]))
    insert!(Ni,1,Ni[1])

    for i in 1:length(e.log_Li)-1
        j=i+1
        if wi_mode=="trapezoidal"
            e.log_wi[j]= logaddexp(e.log_Xi[i], - ((j+1)/Ni[j])) - log(2) #log width of prior mass spanned by the last step-trapezoidal approx
        elseif wi_mode=="simple"
            e.log_wi[j]= logaddexp(e.log_Xi[i], - e.log_Xi[j]/Ni[i]) #log width of prior mass spanned by the last step-simple approx
        else
            throw(ArgumentError("Unsupported wi_mode!"))
        end
        e.log_Liwi[j]=lps(e.log_Li[j],e.log_wi[j]) #log likelihood + log width = increment of evidence spanned by iterate
        e.log_Zi[j]=logaddexp(e.log_Zi[i],e.log_Liwi[j])    #log evidence
        #information- dimensionless quantity
        Hj=lps( 
            (exp(lps(e.log_Liwi[j],-e.log_Zi[j])) * e.log_Li[j]), 
            (exp(lps(e.log_Zi[i],-e.log_Zi[j])) * lps(e.Hi[i],e.log_Zi[i])),
            -e.log_Zi[j])
        Hj === -Inf ? (e.Hi[j]=0.) : (e.Hi[j]=Hj)
    end
end
    
#get a posterior kernel density estimate from the ensemble
function posterior_kde(e::GMC_NS_Ensemble, scale=400., bw=1.)
    !e.sample_posterior && throw(ArgumentError("This ensemble is not set to collect posterior samples!"))

    θ_mat=zeros(length(e.models[1].pos),length(e.models)+length(e.posterior_samples))
    weight_vec=zeros(length(e.models)+length(e.posterior_samples))

    for (n,mrec) in enumerate(e.posterior_samples)
        m=deserialize(mrec.path)
        θ_mat[:,n]=m.θ
        weight_vec[n]=e.log_Liwi[n+1]
    end

    for (n,mrec) in enumerate(e.models)
        m=deserialize(mrec.path)
        p=length(e.posterior_samples)
        θ_mat[:,n+p]=m.θ
        weight_vec[n+p]=m.log_Li + lps(e.log_Xi[end], -length(e.models))
    end

    weight_vec.-=(maximum(weight_vec)-scale)

    return kde!(θ_mat,[bw],exp.(weight_vec))
end
