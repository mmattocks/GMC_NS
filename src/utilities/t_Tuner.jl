abstract type GMC_Tuner end

mutable struct τ_Tuner <: GMC_Tuner
    τ::Float64
    τ_history::Vector{Float64}
    α::Float64
    β::Float64
    a::BitVector
    memory::Integer

    function τ_Tuner(τ1, e::GMC_NS_Ensemble)
        new(τ1,
        zeros(0),
        e.GMC_tune_α,
        e.GMC_tune_β,
        falses(0),
        e.GMC_tune_μ)
    end
end

mutable struct τ_PID <: GMC_Tuner
    τ::Float64
    τ_history::Vector{Float64}
    α::Float64
    β::Float64

    a::BitVector 
    a_history::BitVector
    ℵ::Vector{Float64}
    memory::Integer
    Kp::Float64
    Ki::Float64
    Kd::Float64
    
    ϵ::Float64
    pterm::Float64
    iterm::Float64
    dterm::Float64

    function τ_PID(e::GMC_NS_Ensemble)
        new(e.GMC_init_τ,
        zeros(0),
        e.GMC_tune_α,
        1.,
        falses(0),
        falses(0),
        Vector{Float64}(),
        e.GMC_tune_μ,
        e.GMC_tune_PID...,
        0.,0.,0.,0.)
    end
end

function process_report!(t::τ_PID, report::Bool, update::Bool=true)
    push!(t.a_history, report)
    length(t.a_history)>t.memory ? (t.a=t.a_history[end-t.memory:end]) : (t.a=t.a_history)
    push!(t.ℵ,mean(t.a))
    update && τ_update!(t)
end

function τ_update!(t::τ_PID)
    push!(t.τ_history,t.τ)

    t.ϵ=t.ℵ[end] - t.α
    t.pterm=t.Kp * t.ϵ
    !(t.τ==eps()) && (t.iterm+=(t.Ki * t.ϵ))

    length(t.ℵ) > Int64(round(.1*t.memory)) ? (d = t.ℵ[end-Int64(round(.1*t.memory))] - t.ℵ[end]) : #d is positive if acceptance is declining- in this case want smaller beta and tau
        d=0.
    t.dterm=t.Kd * d

    t.β=max(eps(), 1. + t.pterm + t.iterm - t.dterm)
    t.τ=max(eps(), t.τ * t.β)
end


function init_tune(e::GMC_NS_Ensemble)
    @info "Performing initial GMC timestep tuning on ensemble..."
    τ=1.; tuned=false; τd=1.1; accept=1.; it=1
    while !tuned && it<e.GMC_tune_ι
        test_report=falses(0)
        println(stdout, "τ: $τ, It: $it, Accept: $accept")

        while length(test_report)<Int64(floor(length(e.models)*e.GMC_tune_ρ))
            m=deserialize(rand([m.path for m in e.models]))
            candidate=galilean_trajectory_sample!(m,e,τ)

            candidate.id==m.id ? push!(test_report,0) : push!(test_report,1)
            candidate=0
            GC.gc()
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
        move_cursor_up_while_clearing_lines(stdout, 1)
    end
    println("Done. τ=$τ")
    return τ
end

function tune_τ!(t::τ_Tuner, step_report::BitVector)
    push!(t.τ_history,t.τ)
    t.a=vcat(t.a,step_report)
    ls=length(t.a)
    ls>t.memory && (t.a=t.a[ls-t.memory:end])

    accept=sum(t.a)/ls
    if accept==1.
        t.τ=t.τ*(1+t.β)
    elseif accept<t.α
        t.τ=t.τ*(1-t.β)
    end
end

function Base.show(io::IO, t::τ_PID; progress=false)
    try
        plot=lineplot(t.ℵ, title="Recent Unobstructed Moves", xlabel="Proposals", ylabel="Rate", color=:white, name="Free moves")

        lbls=plot.labels_left
        (arrow="")
        all(t.α .> [v for v in parse.(Float64,values(lbls))]) && (arrow="↑")
        all(t.α .< [v for v in parse.(Float64,values(lbls))]) && (arrow="↓")

        lineplot!(plot, [t.α for i in 1:length(t.ℵ)], name="α setpoint "*arrow, color=:red)
        show(io, plot)
        println()
        printstyled("τ:$(t.τ), α:$(t.α), β:$(t.β)",bold=true)
        println()
        printstyled("Kp:$(t.Kp), Ki:$(t.Ki), Kd:$(t.Kd)",color=:green)
        println()
        printstyled("ptrm:$(t.pterm), itrm:$(t.iterm), dtrm:$(t.dterm)",color=:magenta)
        println()


        progress && return nrows(plot.graphics)+9
    catch
        printstyled("τ_PID NOT AVAILABLE. STANDBY", bold=true)
        progress && return 1
    end
end

function get_tuner_dict(e)
    d=Dict{Int64,τ_PID}()
    for rec in e.models
        d[rec.trajectory]=τ_PID(e)
    end
    return d
end
