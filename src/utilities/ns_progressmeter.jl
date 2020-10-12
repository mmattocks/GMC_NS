const CONVERGENCE_MEMORY=500

mutable struct GMC_NS_Progress{T<:Real} <: AbstractProgress
    interval::T
    dt::AbstractFloat
    start_it::Integer
    counter::Integer
    triggered::Bool
    tfirst::AbstractFloat
    tlast::AbstractFloat
    tstp::AbstractFloat
    printed::Bool        # true if we have issued at least one status update
    desc::AbstractString # prefix to the percentage, e.g.  "Computing..."
    color::Symbol        # default to green
    output::IO           # output stream into which the progress is written
    numprintedvalues::Integer   # num values printed below progress in last iteration
    offset::Integer             # position offset of progress bar (default is 0)

    e::GMC_NS_Ensemble
    tuner::τ_Tuner
    top_m::GMC_NS_Model    

    mean_stp_time::AbstractFloat

    upper_displays::Vector{Vector{Function}}
    lower_displays::Vector{Vector{Function}}
    display_rotate_iterates::Int64
    upperidx::Integer
    loweridx::Integer

    convergence_history::Vector{Float64}
    time_history::Vector{Float64}

    function GMC_NS_Progress{T}(e::GMC_NS_Ensemble,
                               top_m::GMC_NS_Model,
                               tuner::τ_Tuner,
                               interval::T;
                               dt::Real=0.1,
                               desc::AbstractString="Nested Sampling::",
                               color::Symbol=:green,
                               output::IO=stdout,
                               offset::Int=0,
                               start_it::Int=1,
                               upper_displays=[[evidence_display]],
                               lower_displays=[[ensemble_display]],
                               disp_rot_its=0,
                               ) where T
        tfirst = tlast = time()
        printed = false
        new{T}(interval,
         dt,
         start_it,
         start_it,
         false,
         tfirst,
         tlast,
         0., 
         printed, 
         desc,
         color,
         output,
         0,
         offset,
         e,
         tuner,
         top_m,
         0.,
         upper_displays,
         lower_displays,
         disp_rot_its,
         1,
         1,
         zeros(CONVERGENCE_MEMORY),
         zeros(CONVERGENCE_MEMORY))
    end
end

function GMC_NS_Progress(e::GMC_NS_Ensemble, tuner::τ_Tuner, interval::Real; dt::Real=0.1, desc::AbstractString="GMC-NS::", color::Symbol=:green, output::IO=stderr, offset::Integer=0, start_it::Integer=1, upper_displays::AbstractVector{<:AbstractVector{Function}}, lower_displays::AbstractVector{<:AbstractVector{Function}}, disp_rot_its::Integer)
    top_m = deserialize(e.models[findmax([model.log_Li for model in e.models])[2]].path)
    
    return GMC_NS_Progress{typeof(interval)}(e, top_m, tuner, interval, dt=dt, desc=desc, color=color, output=output, offset=offset, start_it=start_it, upper_displays=upper_displays, lower_displays=lower_displays, disp_rot_its=disp_rot_its)
end


function update!(p::GMC_NS_Progress, val, thresh; options...)
    p.counter += 1
    
    p.tstp=time()-p.tlast
    popfirst!(p.time_history)
    push!(p.time_history, p.tstp)

    p.interval = val - thresh
    popfirst!(p.convergence_history)
    push!(p.convergence_history, p.interval)    

    p.top_m.id != basename(p.e.models[findmax([model.log_Li for model in p.e.models])[2]].path) && (p.top_m=deserialize(p.e.models[findmax([model.log_Li for model in p.e.models])[2]].path))

    updateProgress!(p; options...)
end

function updateProgress!(p::GMC_NS_Progress; offset::Integer = p.offset, keep = (offset == 0))
    p.offset = offset
    t = time()
    p.display_rotate_iterates > 0 && p.counter%p.display_rotate_iterates == 0 && rotate_displays(p)

    if p.interval <= 0 && !p.triggered
        p.triggered = true
        if p.printed
            p.triggered = true
            dur = durationstring(t-p.tfirst)
            msg = @sprintf "%s Converged. Time: %s (%d iterations). logZ: %s\n" p.desc dur p.counter p.e.log_Zi[end]
            
            print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
            move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues)
            length(p.upper_displays)>0 ? (upper_lines=display_dash(p, p.upper_displays[p.upperidx])) : (upper_lines=0) 
            printover(p.output, msg, :magenta)
            length(p.lower_displays)>0 ? (lower_lines=display_dash(p, p.lower_displays[p.loweridx])) : (lower_lines=0)
            
            p.numprintedvalues=upper_lines + lower_lines + 1

            if keep
                println(p.output)
            else
                print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))
            end
        end
        return
    end

    if t > p.tlast+p.dt && !p.triggered
        p.counter < CONVERGENCE_MEMORY ? mean_step_time=mean(p.time_history[end-(p.counter-1):end]) :
                                            mean_step_time=mean(p.time_history)
        msg = @sprintf "%s Iterate: %s Recent step time μ: %s Convergence Interval: %g\n" p.desc p.counter hmss(mean_step_time) p.interval

        print(p.output, "\n" ^ (p.offset + p.numprintedvalues))
        move_cursor_up_while_clearing_lines(p.output, p.numprintedvalues)
        length(p.upper_displays)>0 ? (upper_lines=display_dash(p, p.upper_displays[p.upperidx])) : (upper_lines=0)
        printover(p.output, msg, p.color)
        length(p.lower_displays)>0 ? (lower_lines=display_dash(p, p.lower_displays[p.loweridx])) : (lower_lines=0)

        p.numprintedvalues=upper_lines + lower_lines + 1
        print(p.output, "\r\u1b[A" ^ (p.offset + p.numprintedvalues))

        # Compensate for any overhead of printing. This can be
        # especially important if you're running over a slow network
        # connection.
        p.tlast = t + 2*(time()-t)
        p.printed = true
    end
end
                function rotate_displays(p)
                    p.upperidx < length(p.upper_displays) ? p.upperidx+=1 : p.upperidx=1
                    p.loweridx < length(p.lower_displays) ? p.loweridx+=1 : p.loweridx=1
                end

                function hmss(dt)
                    dt<0 ? (dt=-dt; prfx="-") : (prfx="")
                    isnan(dt) && return "NaN"
                    (h,r) = divrem(dt,60*60)
                    (m,r) = divrem(r, 60)
                    (isnan(h)||isnan(m)||isnan(r)) && return "NaN"
                    string(prfx,Int(h),":",Int(m),":",Int(ceil(r)))
                end

                function display_dash(p::GMC_NS_Progress, funcvec::Vector{Function})
                    lines=0
                    for func in funcvec
                        lines+=func(p)
                    end
                    return lines
                end