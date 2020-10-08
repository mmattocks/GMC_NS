"""
    GMC_NS_Ensemble

Abstract supertype for model ensembles intended to be subjected to Galilean Monte Carlo Nested Sampling ('GMC_NS'). Subtyped ensemble structs must implement the following fields as their interface with 'GMC_NS':

    mutable struct Minimal_Ensemble <: GMC_NS_Ensemble
        path::String #ensemble records serialized here
        
        models::Vector{<:GMC_NS_Model_Record} #model paths & lhs
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

        #GMC Sampler settings
        GMC_tune_ι::Int64 #maximum number of initial τ tuning steps
        GMC_tune_ρ::Float64 #apply this scale parameter to the ensemble size to obtain the number of GMC samples per initial tuning step
        GMC_tune_α::Float64 #tune τ such that the acceptance rate is between this value and 1.
        GMC_timestep_η::Float64 #>0. to apply an perturbation term randomly drawn from Normal(0, η*τ)
        GMC_exhaust_σ::Float64 #apply this scale parameter to the ensemble size to obtain the maximum number of GMC iterates to perform before terminating nested sampling
        GMC_reflect_η::Float64 #>0. to apply perturbation to reflection vectors, keep value small

    end 

"""
abstract type GMC_NS_Ensemble end


"""
    show(e<:GMC_NS_Ensemble)
    show(io,e<:GMC_NS_Ensemble)

Display the ensemble path, likelihood histogram, contour, maximum log likelihood, and log evidence (marginal likelihood). 
"""
function Base.show(io::IO, e<:GMC_NS_Ensemble; progress=false)
	livec=[model.log_Li for model in e.models]
	maxLH=maximum(livec)
	printstyled(io, "$(typeof(e)) @ $(e.path)\n", bold=true)
	msg = @sprintf "Contour: %3.6e MaxLH:%3.3e log Evidence:%3.6e" e.contour maxLH e.log_Zi[end]
	println(io, msg)
	hist=UnicodePlots.histogram(livec, title="Ensemble Likelihood Distribution")
	show(io, hist)
	println()
	progress && return(nrows(hist.graphics)+6)
end