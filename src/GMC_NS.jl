module GMC_NS
    using BioBackgroundModels, BioSequences, Distributed, Distributions, Serialization, UnicodePlots
    import ProgressMeter: AbstractProgress, Progress, @showprogress, next!, move_cursor_up_while_clearing_lines, printover, durationstring
    import Printf: @sprintf
    import StatsFuns: logaddexp, logsumexp
    import Random: rand, seed!, shuffle!
    import Distances: euclidean

    include("ensemble/GMC_NS_Ensemble.jl")
    include("model/GMC_NS_Model.jl")
    include("GMC/galilean_trajectory.jl")
    include("GMC/GMC_sample.jl")
    include("nested_sampler/nested_step.jl")
    include("nested_sampler/converge_ensemble.jl")
    include("utilities/ns_progressmeter.jl")
    include("utilities/t_Tuner.jl")
end # module
