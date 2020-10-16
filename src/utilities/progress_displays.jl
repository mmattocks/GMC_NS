function tuning_display(p)
    lines=show(p.output, p.tuner, progress=true)
    println();
    return lines
end

function convergence_display(p)
    ciplot=lineplot([p.counter-(length(p.convergence_history)-1):p.counter...], p.convergence_history, title="Convergence Interval Recent History", xlabel="Iterate",ylabel="CI", color=:yellow)
    lines=nrows(ciplot.graphics)+5
    show(p.output, ciplot); println()
    return lines
end

function evidence_display(p)
    try
        evplot=lineplot(p.e.log_Zi[2:end], title="Evidence History", xlabel="Iterate", color=:red, name="Ensemble logZ")
        lines=nrows(evplot.graphics)+5
        show(p.output, evplot); println()
        return lines
    catch
        printstyled(p.output, "EVIDENCE PLOT UNAVAILABLE. STANDBY\n", bold=true, color=:red); println()
        return 1
    end
end
    
function info_display(p)
    try
        infoplot=lineplot(p.e.Hi[2:end], title="Information History", xlabel="Iterate", color=:green, name="Ensemble H")
        lines=nrows(infoplot.graphics)+5
        show(p.output, infoplot); println()
        return lines
    catch
        printstyled(p.output, "INFORMATION PLOT UNAVAILABLE. STANDBY\n", bold=true, color=:green); println()
        return 1
    end
end

function lh_display(p)
    try
        lhplot=lineplot(p.e.log_Li[2:end], title="Contour History", xlabel="Iterate", color=:magenta, name="Ensemble logLH")
        lines=nrows(lhplot.graphics)+5
        show(p.output, lhplot); println()
        return lines
    catch
        printstyled(p.output, "CONTOUR HISTORY UNAVAILABLE. STANDBY\n", bold=true, color=:magenta); println()
        return 1
    end
end

function liwi_display(p)
    try
        liwiplot=lineplot([max(2,p.counter-(CONVERGENCE_MEMORY-1)):p.counter...],p.e.log_Liwi[max(2,end-(CONVERGENCE_MEMORY-1)):end], title="Recent iterate evidentiary weight", xlabel="Iterate", name="Ensemble log Liwi", color=:cyan)
        lines=nrows(liwiplot.graphics)+5
        show(p.output, liwiplot); println()
        return lines
    catch
        printstyled(p.output, "EVIDENTIARY HISTORY UNAVAILABLE. STANDBY\n", bold=true, color=:cyan); println()
        return 1
    end
end

function ensemble_display(p)
    return lines=show(p.output, p.e, progress=true)
end

function model_display(p)
    println("Current MAP model:")
    return lines=show(p.output, p.top_m, progress=true)
end

function model_obs_display(p)
    println("Current MAP model")
    return lines=show(p.output, p.top_m, p.e, progress=true)
end