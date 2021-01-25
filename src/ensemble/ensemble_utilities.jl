function ensemble_history(e::GMC_NS_Ensemble, bins=25)
    !e.sample_posterior && throw(ArgumentError("This ensemble has no posterior samples to show a history for!"))
    livec=vcat([model.log_Li for model in e.models],[model.log_Li for model in e.posterior_samples])
    show(histogram(livec, nbins=bins))
end

function e_backup(e::GMC_NS_Ensemble, tuner::Dict{Int64,<:GMC_Tuner})
    serialize(string(e.path,'/',"ens"), e)
    serialize(string(e.path,'/',"tuner"), tuner)
end

function clean_ensemble_dir(e::GMC_NS_Ensemble, model_pad::Integer; ignore_warn=false)
    !ignore_warn && e.sample_posterior && throw(ArgumentError("Ensemble is set to retain posterior samples and its directory should not be cleaned!"))
    for file in readdir(e.path)
        !(file in vcat([basename(model.path) for model in e.models],"ens",[string(number) for number in e.model_counter-length(e.models)-model_pad:e.model_counter-1])) && rm(e.path*'/'*file)
    end
end

function complete_evidence(e::GMC_NS_Ensemble)
    return final_logZ = logaddexp(e.log_Zi[end], (logsumexp([model.log_Li for model in e.models] .+  (e.log_Xi[length(e.log_Li)] - log(length(e.models))))))
end

function measure_evidence(e::GMC_NS_Ensemble)
    return ms=measurement(complete_evidence(e),sqrt(abs(e.Hi[end])/length(e.models)))
end


function reset_ensemble!(e::GMC_NS_Ensemble)
    new_e=deepcopy(e)
    for i in 1:length(e.models)
        if string(i) in [basename(record.path) for record in e.models]
            new_e.models[i]=e.models[findfirst(isequal(string(i)), [basename(record.path) for record in e.models])]
        else
            new_e.models[i]=e.posterior_samples[findfirst(isequal(string(i)), [basename(record.path) for record in e.posterior_samples])]
        end
    end

    new_e.contour=minimum([record.log_Li for record in new_e.models])

    new_e.log_Li=[new_e.log_Li[1]]
    new_e.log_Xi=[new_e.log_Xi[1]]
    new_e.log_wi=[new_e.log_wi[1]]
    new_e.log_Liwi=[new_e.log_Liwi[1]]
    new_e.log_Zi=[new_e.log_Zi[1]]
    new_e.Hi=[new_e.Hi[1]]

    new_e.posterior_samples=Vector{typeof(e.models[1])}()

    new_e.model_counter=length(new_e.models)+1

    clean_ensemble_dir(new_e, 0; ignore_warn=true)
    isfile(e.path*"/tuner") && rm(e.path*"/tuner")
    serialize(e.path*"/ens", new_e)

    return new_e
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

function rewind_ensemble(e::GMC_NS_Ensemble,rewind_idx)
    !e.sample_posterior && throw(ArgumentError("An ensemble not retaining posterior samples cannot be rewound!"))
    rewind_idx >= length(e.log_Li) && throw(ArgumentError("rewind_idx must be less than the current iterate!"))

    n=length(e.models)
    max_model_no=length(e.log_Li)+length(e.models)-1
    rewind_model_no=rewind_idx+length(e.models)-1
    new_e = deepcopy(e)

    rm_models=[string(name) for name in rewind_model_no+1:max_model_no]

    filter!(model->!(basename(model.path) in rm_models),new_e.models)
    filter!(model->!(basename(model.path) in rm_models),new_e.posterior_samples)

    while length(new_e.models) < n
        push!(new_e.models,pop!(new_e.posterior_samples))
    end
    new_e.contour=new_e.log_Li[rewind_idx]
    new_e.log_Li=new_e.log_Li[1:rewind_idx]
    new_e.log_Xi=new_e.log_Xi[1:rewind_idx]
    new_e.log_wi=new_e.log_wi[1:rewind_idx]
    new_e.log_Liwi=new_e.log_Liwi[1:rewind_idx]
    new_e.log_Zi=new_e.log_Zi[1:rewind_idx]
    new_e.Hi=new_e.Hi[1:rewind_idx]

    new_e.model_counter=length(new_e.models)+rewind_idx

    return new_e
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
function rectify_ensemble!(e::GMC_NS_Ensemble, wi_mode="trapezoidal")
    Ni=round.(collect(-1:-1:-length(e.log_Li)+1)./e.log_Xi[2:end])
    push!(Ni,Ni[end])

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

function rectify_balance!(o_ens::GMC_NS_Ensemble)
    ens=deepcopy(o_ens)
    println("Rectifying balance for ensemble at $(ens.path)")

    Nvec=Vector{Int64}()
    for i in 2:length(ens.log_Li)
        Xi=ens.log_Xi[i]
        N=round(-(i-1)/Xi)
        push!(Nvec, N)
    end
    insert!(Nvec,1,Nvec[1])

    for chain in 1:ens.t_counter-1
        println("Rectifying $chain ...")
        ch_pfx=ens.path*'/'*string(chain)*'.'
        total_spls=0
        while isfile(ch_pfx*string(total_spls+1))
            total_spls+=1
        end

        lastm=deserialize(ch_pfx*string(1))
        insert_idxs=Vector{Int64}()
        insert_ms=Vector{GMC_NS_Model}()

        for spl in 2:total_spls
            splm=deserialize(ch_pfx*string(spl))

            if splm.θ!=lastm.θ && splm.v!=lastm.v
                push!(insert_idxs,spl)

                m_name=string(chain)*string(spl)

                new_mod=deepcopy(lastm)
                new_mod.v=splm.v
                new_mod.i=spl
                push!(insert_ms, new_mod)

                insert_record=typeof(ens.models[1])(chain, spl, new_mod.pos, ch_pfx*string(new_mod.i),new_mod.log_Li)

                insert_idx=findfirst(i-> >(i,new_mod.log_Li),ens.log_Li)

                insert!(ens.log_Li, insert_idx, ens.log_Li[insert_idx-1])
                insert!(ens.posterior_samples, insert_idx, insert_record)
                insert!(Nvec, insert_idx, Nvec[insert_idx-1])
            end

            lastm=splm
        end

        save_idxs=zeros(Int64, total_spls)

        spl_diff=0
        for idx in 1:total_spls
            idx in insert_idxs && (spl_diff+=1)
            save_idxs[idx]=idx+spl_diff
        end

        post_paths=[rec.path for rec in ens.posterior_samples]
        mod_paths=[rec.path for rec in ens.models]
        modswitch=false

        for idx in total_spls:-1:1
            old_i=ch_pfx*string(idx)

            recidx=findfirst(pth->occursin(old_i,pth), post_paths)
            recidx===nothing && (recidx=findfirst(pth->occursin(old_i,pth), mod_paths); modswitch=true)

            modswitch ? (rec=ens.models[recidx]) : (rec=ens.posterior_samples[recidx])
            rec.i=save_idxs[idx]
            rec.path=ch_pfx*string(idx)

            model=deserialize(old_i)
            model.i=save_idxs[idx]
            serialize(ch_pfx*string(save_idxs[idx]), model)
        end
         
        for (idx,m) in zip(insert_idxs,insert_ms)
            serialize(ch_pfx*string(idx), m)
        end

    end

    ens.log_Xi=	[0] #ie exp(0) = all of the prior is covered
    ens.log_wi=[-Inf] #w0 = 0
    ens.log_Liwi=[-Inf] #Liwi0 = 0
    ens.log_Zi=[-1e300] #Z0 = 0
    ens.Hi=[0] #H0 = 0,
    
    for i in 1:length(ens.log_Li)-1
        j=i+1
        N=Nvec[i]
        push!(ens.log_Xi, -i/N) #log Xi - crude estimate of the iterate's enclosed prior mass
        push!(ens.log_wi, logaddexp(ens.log_Xi[i], - ((j+1)/N)) - log(2)) #log width of prior mass spanned by the last step-trapezoidal approx
        push!(ens.log_Liwi, lps(ens.log_Li[j],ens.log_wi[j])) #log likelihood + log width = increment of evidence spanned by iterate
        push!(ens.log_Zi, logaddexp(ens.log_Zi[i],ens.log_Liwi[j]))    #log evidence
        #information- dimensionless quantity
        Hj=lps( 
            (exp(lps(ens.log_Liwi[j],-ens.log_Zi[j])) * ens.log_Li[j]), 
            (exp(lps(ens.log_Zi[i],-ens.log_Zi[j])) * lps(ens.Hi[i],ens.log_Zi[i])),
            -ens.log_Zi[j])
        Hj === -Inf ? push!(ens.Hi,0.) : push!(ens.Hi, Hj)
    end

    serialize(ens.path*"/ens", ens)

    return ens
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

function ensemble_filter(ens)
    println("Filtering ensemble at $(ens.path)")

    ens_paths=vcat(
        [m.path for m in ens.models],
        [m.path for m in ens.posterior_samples]
    )

    for chain in 1:ens.t_counter-1
        println("Filtering chain $chain...")
        ch_pfx=ens.path*'/'*string(chain)*'.'
        total_spls=0
        while isfile(ch_pfx*string(total_spls+1))
            total_spls+=1
        end

        for spl in 2:total_spls
            splpath=ch_pfx*string(spl)
            if !(splpath in ens_paths)
                rm(splpath)
            end
        end
    end
end