const TrajectoryCourier=Dict{Int64,RemoteChannel}

function deliver_trajectories(results_chan,ts)
    while true
        wait(results_chan)
        result=take!(results_chan)
        put!(ts[result.trajectory],result)
    end
end

function manage_trajectories(nwk, e, jobs_chan)
    min_trajs=[m.trajectory for m in e.models][sortperm([m.log_Li for m in e.models])][1:nwk]
    take!(jobs_chan); put!(jobs_chan,min_trajs)
end

function sample_trajectory(jobs_chan, results_chan)
    while true
        wait(jobs_chan)
        g 9