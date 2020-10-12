mutable struct τ_Tuner
    τ::Float64
    τ_history::Vector{Float64}
    min_acceptance::Float64
    tune_rate::Float64
    successes::BitVector #(memoryxfunc)
    s_memory::Integer


    function τ_Tuner(τ1, e::GMC_NS_Ensemble)
        new(τ1,
        zeros(0),
        e.GMC_tune_α,
        e.GMC_tune_β,
        falses(0),
        e.GMC_tune_μ)
    end
end

function init_tune(e::GMC_NS_Ensemble)
    @info "Performing initial GMC timestep tuning on ensemble..."
    τ=1.; tuned=false; τd=1.1; accept=1.; it=1
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
            last_accept == 1. && (τd*=1+e.GMC_tune_β)
            last_accept < e.GMC_tune_α && (τd=min(1.000001,1-e.GMC_tune_β*τd))
        elseif accept < e.GMC_tune_α
            τ/=τd
            last_accept == 1. && (τd=min(1.000001,1-e.GMC_tune_β*τd))
            last_accept < e.GMC_tune_α && (τd*=1+e.GMC_tune_β)
        else
            return τ
        end

        it+=1
    end
    return τ
end

function tune_τ!(t::τ_Tuner, step_report::BitVector)
    push!(t.τ_history,t.τ)
    t.successes=vcat(t.successes,step_report)
    ls=length(t.successes)
    ls>t.s_memory && (t.successes=t.successes[ls-t.s_memory:end])

    accept=sum(t.successes)/ls
    if accept==1.
        t.τ=t.τ*(1+t.tune_rate)
    elseif accept<t.min_acceptance
        t.τ=t.τ*(1-t.tune_rate)
    end
end

function Base.show(io::IO, t::τ_Tuner; progress=false)
    plot=lineplot(t.τ_history, title="τ History", xlabel="Iterate", color=:white, name="τ")
    show(io, plot)
    println()
    println("τ:$(t.τ), α:$(t.min_acceptance), β:$(t.tune_rate)")

    progress && return nrows(plot.graphics)+7
end