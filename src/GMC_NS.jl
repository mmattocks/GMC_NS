module GMC_NS
    using Distributed, Distributions, Serialization, UnicodePlots
    import ProgressMeter: AbstractProgress, Progress, @showprogress, next!, move_cursor_up_while_clearing_lines, printover, durationstring
    import Printf: @sprintf
    import Random: rand
    import Statistics: unscaled_covzm, det, eigen
    import Serialization: serialize, deserialize
    import LinearAlgebra: dot, normalize
    import BioBackgroundModels: lps
    import StatsFuns: logaddexp, logsumexp
    import Base: show
    import Measurements: measurement
    import NGRefTools: marginals
    import ConjugatePriors: NormalGamma

    #[GMC_τ_death, GMC_init_τ, GMC_tune_μ, GMC_tune_α, GMC_tune_PID, GMC_timestep_η, GMC_reflect_η, GMC_exhaust_σ]

    GMC_DEFAULTS=Vector{Any}([50,1e-6,.5, 50,.8,(.175,1e-5,.1), .001, .01, 10.])
    export GMC_DEFAULTS

    include("ensemble/GMC_NS_Ensemble.jl")
    export GMC_NS_Ensemble
    include("GMC_NS_Model.jl")
    export GMC_NS_Model, GMC_NS_Model_Record
    include("GMC/galilean_trajectory.jl")
    include("utilities/t_Tuner.jl")
    include("utilities/ellipsoid.jl")
    include("nested_sampler/search_patterns.jl")
    include("nested_sampler/nested_step.jl")
    include("nested_sampler/converge_ensemble.jl")
    export converge_ensemble!
    include("utilities/coordinate_utils.jl")
    export box_bound!, to_prior, to_unit_ball, box_reflect!
    include("utilities/ensemble_utilities.jl")
    export ensemble_history, clean_ensemble_dir, complete_evidence, reset_ensemble!, move_ensemble!, copy_ensemble!, rewind_ensemble, show_models, show_models_e, get_model
    include("utilities/ns_progressmeter.jl")
    include("utilities/progress_displays.jl")
    export tuning_display,evidence_display,convergence_display,info_display,lh_display,liwi_display,ensemble_display,model_display,model_obs_display
    include("normal/Normal_Model.jl")
    include("normal/Normal_Ensemble.jl")
    export Normal_Ensemble
    include("normal/LogNormal_Model.jl")
    include("normal/LogNormal_Ensemble.jl")
    export LogNormal_Ensemble
end # module
