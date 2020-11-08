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
        if t.τ > e.GMC_τ_death
            model_selected=galilean_search!(m,e,t,least_likely_idx,samples,sample_limit)
        elseif (N-1)<e.GMC_Nmin #if removing the particle would put the sample under the minimum size, resample
            model_selected=diffusion_search(e,tuners,least_likely_idx,N,e.GMC_τ_death)
            samples+=1
        else #if not, push to the posterior samples and let it die
            model_selected=true
            e.sample_posterior && push!(e.posterior_samples, e.models[least_likely_idx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
            deleteat!(e.models, least_likely_idx) #remove least likely model record
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