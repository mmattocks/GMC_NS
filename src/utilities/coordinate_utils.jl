function to_prior(pos,prior)
    return quantile(prior, .5+.5*pos)
end

function to_unit_ball(pos,prior)
    return (cdf(prior,pos)-.5)/.5
end


function box_bound!(pos,box)
    if !(all(box[:,1].<pos.<box[:,2]))
        pos[pos.<box[:,1]].=box[:,1][pos.<box[:,1]]
        pos[pos.>box[:,2]].=box[:,2][pos.>box[:,2]]
    end
end