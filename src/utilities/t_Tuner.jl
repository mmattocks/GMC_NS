mutable struct τ_Tuner
    successes::BitVector #(memoryxfunc)
    τ::Float64
    τ_history::Vector{Float64}
    τ_Tuner(τ1)=new(falses(0),τ1,zeros(0))
end

function init_tune(e<:GMC_NS_Ensemble)
    @info "Performing initial GMC timestep tuning on ensemble..."
    τ=1.; tuned=false; τd=1.1; accept=1.; it+1
    while !tuned && it<e.GMC_tune_ι
        test_report=falses(0)
        while length(test_report)<Int64(floor(length(e.models)*e.GMC_tune_ρ))
            m=deserialize(rand([m.path for m in e.models]))
            candidate=galilean_trajectory_sample!(m,e,τ)

            candidate.id==m.id ? push!(test_report,0) : push!(test_report,1)
        end
        
        last_accept=accept
        accept=sum(test_report)/length(test_report)

        if accept == 1.
            τ*=τd
            last_accept == 1. && τd*=τd
            last_accept < e.GMC_τune_α && τd=min(1.000001,.9*τd)
        elseif accept < e.GMC_τune_α
            τ/=τd
            last_accept == 1. && τd=min(1.000001,.9*τd)
            last_accept < e.GMC_τune_α && τd*=τd
        else
            return τ
        end

        it+=1
    end
    return t
end

function Base.show(io::IO, tuner::Permute_Tuner; progress=false)
    show(io, tuner.tabular_display, rowlabel=:I, summary=false)
    progress && return(size(tuner.tabular_display,1)+4)
end