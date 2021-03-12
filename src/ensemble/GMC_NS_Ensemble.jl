"""
    GMC_NS_Ensemble

Abstract supertype for model ensembles intended to be subjected to Galilean Monte Carlo Nested Sampling ('GMC_NS'). Subtyped ensemble structs must implement the following fields as their interface with 'GMC_NS':

    mutable struct Minimal_Ensemble <: GMC_NS_Ensemble
        path::String #ensemble records serialized here
        
        models::Vector{<:GMC_NS_Model_Record} #records of model paths & lhs
        model_initλ::Function #GMC_NS_Model constructor

        contour<:AbstractFloat #current log_Li[end], init minimum([record.log_Li for record in models])
        log_Li::Vector{<:AbstractFloat} #likelihood of lowest-ranked model at iterate i, init [-Inf]
        log_Xi::Vector{<:AbstractFloat} #amt of prior mass included in ensemble contour at Li, init [0]
        log_wi::Vector{<:AbstractFloat} #width of prior mass covered in this iterate, init [-Inf]
        log_Liwi::Vector{<:AbstractFloat} #evidentiary weight of the iterate, init [-Inf]
        log_Zi::Vector{<:AbstractFloat} #ensemble evidence, init [typemin(FloatType)]
        Hi::Vector{<:AbstractFloat} #ensemble information, init [0]

        obs::Any #observations formatted for GMC_NS_Model likelihood Function
        priors::Vector{<:Distribution} #prior distributions on GMC_NS_Model parameters, formatted for model_initλ
        constants::Vector #fixed parameters to be passed to GMC_NS_Model constructor

        sample_posterior::Bool #switch for posterior_samples collection
        posterior_samples::Vector{<:GMC_NS_Model_Record}

        #GMC tuning settings
        GMC_tune_ι<:Integer #maximum number of initial τ tuning steps
        GMC_tune_ρ<:AbstractFloat #use ensemble size*ρ GMC samples per initial tuning step
        GMC_tune_μ<:Integer #tune τ on the basis of this number of the most recent steps
        GMC_tune_α<:AbstractFloat #tune τ such that the acceptance rate is between this value and 1.
        GMC_tune_β<:AbstractFloat #multiply τ by 1+β (sampling acceptance rate=1.) or 1-β (rate <α)

        #GMC Sampler settings
        GMC_timestep_η<:AbstractFloat #>0. to apply an perturbation term randomly drawn from Normal(0, η*τ) to τ at each step (can assist equilibration)
        GMC_reflect_η<:AbstractFloat #>0. to apply perturbation to reflection vectors, keep value small (can assist equilibration)
        GMC_exhaust_σ<:AbstractFloat #apply this scale parameter to the ensemble size to obtain the maximum number of GMC iterates to perform before terminating nested sampling
        GMC_chain_κ<:Integer #maximum chain length

        t_counter<:Integer #counter for next instantiated trajectory number
    end 

"""
abstract type GMC_NS_Ensemble end


"""
    show(e<:GMC_NS_Ensemble)
    show(io,e<:GMC_NS_Ensemble)

Display the ensemble path, likelihood histogram, contour, maximum log likelihood, and log evidence (marginal likelihood). 
"""
function Base.show(io::IO, e::GMC_NS_Ensemble; progress=false)
	livec=[model.log_Li for model in e.models]
	maxLH=maximum(livec)
	printstyled(io, "$(typeof(e)) @ $(e.path)\n", bold=true)
	msg = @sprintf "Contour: %3.6e MaxLH:%3.3e log Evidence:%3.6e" e.contour maxLH e.log_Zi[end]
    println(io, msg)
    if length(unique(livec)) > 1
        try
            hist=UnicodePlots.histogram(livec, title="Ensemble Likelihood Distribution")
            show(io, hist)
            println()
            progress && return(nrows(hist.graphics)+6)
        catch
            printstyled("HISTOGRAM UNAVAILABLE", color=:green)
            println()
            progress && return 3
        end
    else
        printstyled("HISTOGRAM UNAVAILABLE", color=:green)
        println()
        progress && return 3
    end
end