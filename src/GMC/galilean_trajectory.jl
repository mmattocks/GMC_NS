"""
    galilean_trajectory_sample(e, m, τ)

Given 'm <: GMC_NS_Model', 'e <: GMC_NS_Ensemble', and time step τ, attempt to find new model obeying 'new_m.log_Li > e.contour', first by permuting 'm.θ' by the distance given by velocity 'm.v' and time step 'τ', with 'm' retaining particle velocity v:

''[θ,v] -> [θ′=θ+τv, v]''

If 'e.GMC_timestep_η' is a positive value, τ has a random error term of Normal(0,τ*η) applied.

 If this fails, the rejected particle "reflects" off the contour boundary as represented by the gradient normal 'n':

 ''[θ,v] -> [θ′′=θ+τv+τv′, v′=v-2n(n∘v)]''

If 'e.GMC_reflect_η' is a positive value, this value (η below) is used with a newly sampled isotropic vector to calculate a perturbed reflection vector 'v′′':

''v′′=v′cosη+rsinη, r=Normal(0,I)''
 
Finally, if reflection fails to produce a new particle inside the contour, invert the particle's velocity and return it to the ensemble, preserving detailed balance:

''[θ,v] -> [θ, v′=-v]''

"""
function galilean_trajectory_sample!(m, e, tuner, cache)
    t=m.trajectory; next_i=m.i+1; cache_i=m.i+2

    e.GMC_timestep_η > 0. ? (τ=max(tuner.τ+rand(Normal(0,tuner.τ*e.GMC_timestep_η)),eps())) : (τ=tuner.τ) #perturb τ if specified
    d=τ*m.v #distance vector for given timestep
    fwd_pos,adj_d,box_rflct=box_move(m.pos,d,e.box)
    
    if !box_rflct #if we're not stopped by the sampling box
        fwd_m=e.model_initλ(t, next_i, to_prior.(fwd_pos,e.priors), fwd_pos, m.v, e.obs, e.constants...)  #try to proceed along distance vector
        if m.θ==fwd_m.θ #if forward motion is no longer making a difference to the parameter vector (ie timestep is too small), kill the trajectory and bail out
            tuner.τ=e.GMC_τ_death
            return m, nothing
        end
    else
        fwd_m=m #placeholder, no effect on box reflection
    end

    if box_rflct || (fwd_m.log_Li < m.log_Li) #if no model is available from the position along the distance vector, search for reflections and store in cache
        process_report!(tuner, false)

        lhδ=lps(m.log_Li,-fwd_m.log_Li) #get the likelihood difference between the two points
        n=boundary_norm(adj_d,lhδ) #get the gradient normal for the distance and lhδ
        box_rflct ? (v′=box_reflect(m.pos,e.box,m.v)) : v′=reflect(m.v,n,e.GMC_reflect_η) #get the reflected velocity vector off the boundary "east", or off the sampling box
        rd=τ*v′ #reflected distance given the timestep
        east_pos,_=box_move(m.pos,rd,e.box)
        cache=e.model_initλ(t, cache_i, to_prior.(east_pos,e.priors), east_pos, v′, e.obs, e.constants...) #try going east, or reflecting off the box
        if cache.log_Li < m.log_Li && !box_rflct #no west reflection off the sampling box
            process_report!(tuner, false)

            west_pos,_=box_move(m.pos,-rd,e.box)
            cache=e.model_initλ(t, cache_i, to_prior.(west_pos,e.priors), west_pos, -v′, e.obs, e.constants...)  #if east reflection fails, try west
        end
        if cache.log_Li < m.log_Li
            process_report!(tuner, false)
            south_pos,_=box_move(m.pos,-d,e.box)
            cache=e.model_initλ(t, cache_i, to_prior.(south_pos,e.priors),south_pos,-m.v, e.obs, e.constants...)  #if west or box reflection fails, try south (reversed along the particle's vector)
        end
        cache.log_Li < m.log_Li && (process_report!(tuner, false)) #if that fails wait for smaller τ
    end

    if cache!==nothing && cache.log_Li >= m.log_Li #a reflection that results in a model above contour exists, so return a model in place with the appropriate velocity vector and the cached model down the reflection vector
        r_m=construct_reflection(m, next_i, cache.v)
        return (r_m, cache)
    else #if no reflections occurred or no reflected positions above the contour were found - either forward search succeeded or the whole search failed
        cache!==nothing && (cache=nothing) #reset cache if it's been set with a model below contour
        fwd_m.i == next_i && fwd_m.log_Li >= m.log_Li ? (process_report!(tuner, true); return fwd_m, cache) : (return m, cache) #if the forward model is a new iterate and at or above the contour, return that, otherwise return the original model
    end
end

                function box_move(pos,d,box)
                    if !(all(box[:,1].<pos+d.<box[:,2])) #if not all coordinates after move would be in the box
                        boundary_d=box.-pos #get distances to box boundaries
                        boundary_t=boundary_d./d #get time to hit boundaries, given particle displacment d

                        for idx in findall(x->x>=0.,boundary_t) #for any boundaries the particle is riding (distance 0), give a small positive time if the particle is headed into the boundary to allow that dimension to be selected for bounding
                            (dim,bound)=(idx[1],idx[2])
                            boundary_t[idx]==0. && ((bound==1 && d[dim] < 0.) || (bound==2 && d[dim] > 0.)) && (boundary_t[idx]=eps())
                        end

                        first_b_idx=findfirst(isequal(minimum(boundary_t[findall(x->x>0.,boundary_t)])),boundary_t) #find the index of the first boundary that would be hit
                        first_b_d=boundary_d[first_b_idx]
                        d=d.*(first_b_d/d[first_b_idx[1]])
                    end

                    all(d.==0.) ? box_reflect=true : box_reflect=false

                    return pos+d, d, box_reflect
                end

                function box_reflect(pos,box,v)
                    b=isapprox.(pos,box)
                    low_idxs=[v[i]<0&&b[i,1] for i in 1:length(v)]
                    hi_idxs=[v[i]>0&&b[i,2] for i in 1:length(v)]

                    rv=copy(v)
                    rv[low_idxs]=-v[low_idxs]
                    rv[hi_idxs]=-v[hi_idxs]
                    return rv
                end

                function boundary_norm(v, lhδ)
                    b=lhδ./v
                    b[b.==Inf].=prevfloat(Inf) #rectify overflow
                    b[b.==0.].=rand([nextfloat(0.),prevfloat(0.)]) #rectify exact 0.
                    b[b.==-Inf].=nextfloat(-Inf) #rectify underflow
                    return normalize(b)
                end

                function reflect(v, n, η)
                    v′=v-(2*n*dot(n,v))
                    η > 0. ? (return r_perturb(v′,η)) : (return v′)
                end

                function r_perturb(v′, η)
                    return v′*cos(η) + rand(MvNormal(length(v′),1.))*sin(η)
                end