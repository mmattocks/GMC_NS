function converge_ensemble!(e::GMC_NS_Ensemble; max_iterates=typemax(Int64), τ0::Float64=0., backup::Tuple{Bool,Integer}=(false,0), clean::Tuple{Bool,Integer,Integer}=(false,0,0), verbose::Bool=false, converge_criterion::String="standard", converge_factor::AbstractFloat=.01, τPID::NTuple{3,Float64}=(.1,1e-5,.1), progargs...)
    N = length(e.models); curr_it=length(e.log_Li)

    if curr_it==1 || !isfile(e.path*"/tuner")
        τ0 > 0. ? (init_τ=τ0) : (init_τ=init_tune(e))
        tuner=τ_PID(init_τ, e, τPID...);
    else
        tuner=deserialize(e.path*"/tuner") #restore tuner from saved if any
    end

    meter = GMC_NS_Progress(e, tuner, 0.; start_it=curr_it, progargs...)

    converge_check = get_convfunc(converge_criterion)
    while !converge_check(e, converge_factor) && (curr_it <= max_iterates)
        warn = nested_step!(e, tuner)
        warn == 1 && (@error "Failed to find new models, aborting at current iterate."; return e)
        curr_it += 1

        backup[1] && curr_it%backup[2] == 0 && e_backup(e,tuner) #every backup interval, serialise the ensemble and instruction
        clean[1] &&  !e.sample_posterior && curr_it%clean[2] == 0 && clean_ensemble_dir(e,clean[3]) #every clean interval, remove old discarded models

        update!(meter, converge_check(e,converge_factor,vals=true)...)
    end

    if converge_check(e,converge_factor)
        final_logZ = complete_evidence(e)
        @info "Job done, sampled to convergence. Final logZ $final_logZ"

        e_backup(e,tuner)
        clean[1] && !e.sample_posterior && clean_ensemble_dir(e,0) #final clean
        return final_logZ
    elseif curr_it==max_iterates
        @info "Job done, sampled to maximum iterate $max_iterates. Convergence criterion not obtained."

        e_backup(e,tuner)
        clean[1] && !e.sample_posterior && clean_ensemble_dir(e,0) #final clean
        return e.log_Zi[end]
    end
end

function converge_ensemble!(e::GMC_NS_Ensemble, wk_pool::Vector{Int64}; max_iterates=typemax(Int64), τ0::Float64=0., backup::Tuple{Bool,Integer}=(false,0), clean::Tuple{Bool,Integer,Integer}=(false,0,0), verbose::Bool=false, converge_criterion::String="standard", converge_factor::AbstractFloat=.001, progargs...)
    N = length(e.models)
    converge_check = get_convfunc(converge_criterion)

    curr_it=length(e.log_Li)
    if curr_it==1 || !isfile(e.path*"/tuner")
        τ0 > 0. ? (init_τ=τ0) : (init_τ=init_tune(e))
        tuner=τ_Tuner(init_τ, e);
    else
        tuner=deserialize(e.path*"/tuner") #restore tuner from saved if any
    end

    @info "wrreeowow"

    model_chan= RemoteChannel(()->Channel{GMC_NS_Model}(Inf)) #channel to take EM iterates off of
    job_chan = RemoteChannel(()->Channel{Tuple{<:AbstractVector{<:GMC_NS_Model_Record}, Float64, Float64}}(1))
    put!(job_chan,(e.models, e.contour, tuner.τ))    
    if !converge_check(e,converge_factor) #sequence workers only if not already converged
        for wk in wk_pool
            @info "wk"
            remote_do(println, wk, galilean_trajectory_sample!(e,job_chan,model_chan))
            @info "jojojo"
            remote_do(galilean_trajectory_sample!,wk,e,job_chan,model_chan)
            @info "hehehehe"
        end
        #@async sequence_workers(wk_pool, galilean_trajectory_sample!, e, job_chan, model_chan)
    end

    meter = GMC_NS_Progress(e, tuner, 0.; start_it=curr_it, progargs...)

    while !converge_check(e, converge_factor) && (curr_it <= max_iterates)
        warn, step_report = nested_step!(e, model_chan) #step the ensemble
        warn == 1 &&
                (@error "All workers failed to find new models, aborting at current iterate."; return e) #if there is a warning, iust return the ensemble and print info
        curr_it += 1        
        
        tune_τ!(tuner, step_report)

        backup[1] && curr_it%backup[2] == 0 && e_backup(e,tuner) #every backup interval, serialise the ensemble and instruction
        clean[1] &&  !e.sample_posterior && curr_it%clean[2] == 0 && clean_ensemble_dir(e,clean[3]) #every clean interval, remove old discarded models

        update!(meter, converge_check(e,converge_factor,vals=true)...)
    end

    take!(job_chan); put!(job_chan, (e.models, e.contour, 0.)) #0. τ terminates worker functions

    if converge_check(e,converge_factor)
        final_logZ = complete_evidence(e)
        @info "Job done, sampled to convergence. Final logZ $final_logZ"

        e_backup(e,tuner)
        clean[1] && !e.sample_posterior && clean_ensemble_dir(e,0) #final clean
        return final_logZ
    elseif curr_it==max_iterates
        @info "Job done, sampled to maximum iterate $max_iterates. Convergence criterion not obtained."

        e_backup(e,tuner)
        clean[1] && !e.sample_posterior && clean_ensemble_dir(e,0) #final clean
        return e.log_Zi[end]
    end
end

                function evidence_converge(e, evidence_fraction; vals=false)
                    val=lps(findmax([model.log_Li for model in e.models])[1],  e.log_Xi[end])
                    thresh=lps(log(evidence_fraction),e.log_Zi[end])

                    isinf(val)||isinf(thresh) && return false

                    vals ? (return val, thresh) : (return val<thresh)
                end

                function compress_converge(e, compression_ratio; vals=false)
                    val=findmax([model.log_Li for model in e.models])[1]-e.contour
                    thresh=compression_ratio

                    isinf(val) && return false
                    vals ? (return val, thresh) : (return val<thresh)
                end

                function get_convfunc(criterion)
                    if criterion == "standard"
                        return evidence_converge
                    elseif criterion == "compression"
                        return compress_converge
                    else
                        throw(ArgumentError("Convergence criterion $criterion not supported! Try \"standard\" or \"compression\"."))
                    end
                end
