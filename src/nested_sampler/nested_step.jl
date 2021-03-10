"""
    nested_step!(e<:GMC_NS_Ensemble)

Step and update ensemble e by sampling the e.contour-bounded prior mass via Galilean Monte Carlo.

"""

function nested_step!(e::GMC_NS_Ensemble, tuners::Dict{Int64,τ_PID}, cache::Union{GMC_NS_Model,Nothing})
    N = length(e.models) #number of sample models/particles on the posterior surface
    i = length(e.log_Li) #iterate number, index for last values
    j = i+1 #index for newly pushed values

    e.contour, least_likely_idx = findmin([model.log_Li for model in e.models])
    Li_model = e.models[least_likely_idx]
    m = deserialize(Li_model.path)
    t=tuners[Li_model.trajectory]

    #SELECT NEW MODEL, SAVE TO ENSEMBLE DIRECTORY, CREATE RECORD AND PUSH TO ENSEMBLE
    
    cache===nothing ? (model_selected=false; cache_insert=false) : (model_selected=true; cache_insert=true)
    samples=1; sample_limit=Int64(floor(length(e.models)*e.GMC_exhaust_σ))

    while !model_selected && samples<=sample_limit
        if t.τ > e.GMC_τ_death && m.i <= e.GMC_chain_κ #while the particle's timestep is greater than the minimum to keep a chain alive, and the chain is below the maximum allowed length, sample by advancing it along a galilean trajectory
            model_selected,cache = galilean_search!(m,e,t,least_likely_idx,samples,sample_limit,cache) #ensemble is modified by galilean_search! directly
        elseif (N-1)<e.GMC_Nmin #if removing the particle would put the sample under the minimum size, resample by diffusion from an existing particle
            model_selected=diffusion_search!(e,tuners,least_likely_idx,N,e.GMC_τ_death) #ensemble modified by diffusion_search!
            samples+=1
        else #if not, push to the posterior samples and let the trajectory chain die
            model_selected=true
            e.sample_posterior && push!(e.posterior_samples, e.models[least_likely_idx])#if sampling posterior, push the model record to the ensemble's posterior samples vector
            deleteat!(e.models, least_likely_idx) #remove least likely model record from the active particles in the ensemble
        end
    end

    if cache_insert #in this case a cached precalculated reflection model can be used
        insert_reflection!(e, cache, least_likely_idx) #insert the model into the ensemble
        process_report!(t, true)
        cache=nothing #reset the cache
    end
        
    samples==sample_limit && (return 1, cache) #failed to find a new sample after sampling up to the limit, return an error code
      
    #UPDATE ENSEMBLE QUANTITIES   
    push!(e.log_Li, minimum([model.log_Li for model in e.models])) #log likelihood of the least likely model - the current ensemble ll contour at Xi
    push!(e.log_Xi, -i/N) #log Xi - crude estimate of the iterate's enclosed prior mass
    push!(e.log_wi, logaddexp(e.log_Xi[i], - ((j+1)/N)) - log(2)) #log width of prior mass spanned by the last step-trapezoidal approx
    push!(e.log_Liwi, lps(e.log_Li[j],e.log_wi[j])) #log likelihood + log width = increment of evidence spanned by iterate
    push!(e.log_Zi, logaddexp(e.log_Zi[i],e.log_Liwi[j]))    #log evidence
    #information- dimensionless quantity. cf. MultiNest @ https://github.com/farhanferoz/MultiNest/blob/master/MultiNest_v3.12/nested.F90- MultiNest now gives info only at the final step
    Hj=lps( 
        (exp(lps(e.log_Liwi[j],-e.log_Zi[j])) * e.log_Li[j]), #information contribution of this step
        (exp(lps(e.log_Zi[i],-e.log_Zi[j])) * lps(e.Hi[i],e.log_Zi[i])), #rescale last information
        -e.log_Zi[j]) 
    Hj === -Inf ? push!(e.Hi,0.) : push!(e.Hi, Hj) #prevent problems from early strings of models with -Inf log likelihoods

    return 0, cache
end