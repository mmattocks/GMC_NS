module GMC_NS
    using  Distributions, Serialization, UnicodePlots
    import ProgressMeter: AbstractProgress, Progress, @showprogress, next!, move_cursor_up_while_clearing_lines, printover, durationstring
    import Printf: @sprintf
    import Random: rand
    import Statistics: unscaled_covzm, det, eigen
    import Serialization: serialize, deserialize
    import LinearAlgebra: dot, normalize
    import BioBackgroundModels: lps
    import StatsFuns: logaddexp, logsumexp, logsubexp
    import Base: show
    import Measurements: measurement
    import ConjugatePriors: NormalGamma, NormalInverseGamma
    import NGRefTools: MarginalTDist
    import KernelDensity: kde, AbstractKDE

    #GMC settings vector fields: [GMC_Nmin::Int64, GMC_τ_death::Float64, GMC_init_τ::Float64, GMC_tune_μ::Int64, GMC_tune_α::Float64, GMC_tune_PID::NTuple{3,Float64}, GMC_timestep_η::Float64, GMC_reflect_η::Float64, GMC_exhaust_σ::Float64, GMC_chain_κ::Int64]

    #should GMC_reflect_η be smaller?
    GMC_DEFAULTS=Vector{Any}([50,1e-5,1., 4,.25,(.3,0.,.0), .25, .005, 10., typemax(Int64)])
    export GMC_DEFAULTS

    include("ensemble/GMC_NS_Ensemble.jl")
    export GMC_NS_Ensemble
    include("model/GMC_NS_Model.jl")
    export GMC_NS_Model, GMC_NS_Model_Record
    include("model/model_utilities.jl")
    include("GMC/galilean_trajectory.jl")
    include("utilities/t_Tuner.jl")
    include("utilities/ellipsoid.jl")
    include("nested_sampler/search_patterns.jl")
    include("nested_sampler/nested_step.jl")
    include("nested_sampler/converge_ensemble.jl")
    export converge_ensemble!
    include("utilities/coordinate_utils.jl")
    export box_bound!, to_prior, to_unit_ball, box_reflect!
    include("ensemble/ensemble_utilities.jl")
    export ensemble_history, clean_ensemble_dir, measure_evidence, reset_ensemble!, move_ensemble!, copy_ensemble, rewind_ensemble, show_models, show_models_e, get_model, rectify_ensemble!, posterior_kde
    include("utilities/ns_progressmeter.jl")
    include("utilities/progress_displays.jl")
    export tuning_display,evidence_display,convergence_display,info_display,lh_display,liwi_display,ensemble_display,model_display,model_obs_display
    include("utilities/stats.jl")
    include("model/normal/Normal_Model.jl")
    include("model/normal/Normal_Ensemble.jl")
    export Normal_Ensemble
    include("model/normal/LogNormal_Model.jl")
    include("model/normal/LogNormal_Ensemble.jl")
    export LogNormal_Ensemble
    include("model/eggbox/Eggbox_Model.jl")
    include("model/eggbox/Eggbox_Ensemble.jl")
    export Eggbox_Ensemble
end # module
