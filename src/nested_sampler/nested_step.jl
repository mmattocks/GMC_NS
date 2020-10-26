"""
    nested_step!(e<:GMC_NS_Ensemble)

Step and update ensemble e by sampling the e.contour-bounded prior mass via Galilean Monte Carlo.

"""

function nested_step!(e::GMC_NS_Ensemble, tuners::Dict{Int64,τ_PID})
    N = length(e.models) #number of sample models/particles on the posterior surface
    i = length(e.log_Li) #iterate number, index for last values
    j = i+1 #index for newly pushed values

    e.contour, least_likely_idx = findmin([model.log_Li for model in e.models])
    Li_model = e.models[least_likely_idx]
    m = deserialize(Li_model.path)
    t=tuners[Li_model.trajectory]

    #SELECT NEW MODEL, SAVE TO ENSEMBLE DIRECTORY, CREATE RECORD AND PUSH TO ENSEMBLE
    model_selected=false
    samples=1; sample_limit=Int64(floor(length(e.models)*e.GMC_exhaust_σ))

    while !model_selected && samples<=sample_limit
        t.τ > e.GMC_τ_death ? (resample=false) : (resample=true)
        if !resample
            model_selected=galilean_search!(m,e,t,least_likely_idx,samples,sample_limit)
        else
            model_selected=ellipsoid_resample!(e,tuners,least_likely_idx,samples,sample_limit)
        end
    end

    samples==sample_limit && (return 1)
      
    #UPDATE ENSEMBLE QUANTITIES   
    push!(e.log_Li, minimum([model.log_Li for model in e.models])) #log likelihood of the least likely model - the current ensemble ll contour at Xi
    push!(e.log_Xi, -i/N) #log Xi - crude estimate of the iterate's enclosed prior mass
    push!(e.log_wi, lps(e.log_Xi[i], -((j+1)/N)-log(2))) #log width of prior mass spanned by the last step-trapezoidal approx
    push!(e.log_Liwi, lps(e.log_Li[j],e.log_wi[j])) #log likelihood + log width = increment of evidence spanned by iterate
    push!(e.log_Zi, logaddexp(e.log_Zi[i],e.log_Liwi[j]))    #log evidence
    #information- dimensionless quantity
    push!(e.Hi, lps(
            (exp(lps(e.log_Liwi[j],-e.log_Zi[j])) * e.log_Li[j]), #term1
            (exp(lps(e.log_Zi[i],-e.log_Zi[j])) * lps(e.Hi[i],e.log_Zi[i])), #term2
            -e.log_Zi[j])) #term3

    return 0
end

function nested_step!(e::GMC_NS_Ensemble, model_chan::RemoteChannel)
    N = length(e.models)+1 #number of sample models/particles on the posterior surface- +1 as one has been removed in the distributed dispatch for converge_ensemble
    i = length(e.log_Li) #iterate number, index for last values
    j = i+1 #index for newly pushed values

    #SELECT NEW MODEL, SAVE TO ENSEMBLE DIRECTORY, CREATE RECORD AND PUSH TO ENSEMBLE
    model_selected=false; report=falses(0);
    samples=0; sample_limit=Int64(floor(length(e.models)*e.GMC_exhaust_σ))

    while !model_selected && samples<=sample_limit
        @async wait(model_chan)
        candidate = take!(model_chan)

        if candidate.id!=0 #particle reversed rather than new sample
            push!(report, false)
            path=e.path*'/'*string(candidate.id)
            serialize(path,candidate)
            samples+=1
        else
            push!(report, true)

            model_selected=true
            candidate.id=e.model_counter
            
            deleteat!(e.models, least_likely_idx) #remove least likely model record
            e.sample_posterior && push!(e.posterior_samples, Li_model)#if sampling posterior, push the model record to the ensemble's posterior samples vector

            new_model_record = Model_Record(string(e.path,'/',e.model_counter), candidate.log_Li)
            push!(e.models, new_model_record) #insert new model record
            serialize(new_model_record.path, candidate)

            e.model_counter +=1
        end
    end

    samples==sample_limit && (return 1, report)
    
    #UPDATE ENSEMBLE QUANTITIES   
    push!(e.log_Li, minimum([model.log_Li for model in e.models])) #log likelihood of the least likely model - the current ensemble ll contour at Xi
    push!(e.log_Xi, -i/N) #log Xi - crude estimate of the iterate's enclosed prior mass
    push!(e.log_wi, lps(e.log_Xi[i], -((j+1)/N)-log(2))) #log width of prior mass spanned by the last step-trapezoidal approx
    push!(e.log_Liwi, lps(e.log_Li[j],e.log_wi[j])) #log likelihood + log width = increment of evidence spanned by iterate
    push!(e.log_Zi, logaddexp(e.log_Zi[i],e.log_Liwi[j]))    #log evidence

    #information- dimensionless quantity
    push!(e.Hi, lps(
            (exp(lps(e.log_Liwi[j],-e.log_Zi[j])) * e.log_Li[j]), #term1
            (exp(lps(e.log_Zi[i],-e.log_Zi[j])) * lps(e.Hi[i],e.log_Zi[i])), #term2
            -e.log_Zi[j])) #term3

    return 0, step_report
end