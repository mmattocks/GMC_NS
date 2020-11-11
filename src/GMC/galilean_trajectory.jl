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
function galilean_trajectory_sample!(m, e, tuner)
    t=m.trajectory;i=m.i+1

    e.GMC_timestep_η > 0. ? (τ=max(tuner.τ+rand(Normal(0,tuner.τ*e.GMC_timestep_η)),eps())) : (τ=tuner.τ) #perturb τ if specified
    d=τ*m.v #distance vector for given timestep
    box_rflct=false
    fwd_pos,adj_d=box_move(m.pos,d,e.box)
    isapprox(fwd_pos,m.pos)&&(box_rflct=true)

    # println(e.box)
    # println(m.pos)
    # println(d)
    # println(fwd_pos)
    
    !box_rflct ? (new_m=e.model_initλ(t, i, to_prior.(fwd_pos,e.priors), fwd_pos, m.v, e.obs, e.constants...)) : (new_m=m) #try to proceed along distance vector

    # println("fwd: $(new_m.log_Li-m.log_Li) ")

    if box_rflct || (new_m.log_Li < m.log_Li) #if that didnt work
        process_report!(tuner, false)

        lhδ=lps(m.log_Li,-new_m.log_Li) #get the lhδ between the two points
        n=boundary_norm(adj_d,lhδ) #get the gradient normal for the distance and lhδ
        box_rflct ? (v′=box_reflect(m.pos,e.box,m.v)) : v′=reflect(m.v,n,e.GMC_reflect_η) #get the reflected velocity vector off the boundary "east"
        rd=τ*v′ #reflected distance given the timestep
        # println("v': $v′, rd $rd")
        east_pos,_=box_move(m.pos,rd,e.box)
        new_m=e.model_initλ(t, i, to_prior.(east_pos,e.priors), east_pos, v′, e.obs, e.constants...) #try going east
        # println("newθ: $(to_prior.(east_pos,e.priors))")

        # println("east: $(new_m.log_Li-m.log_Li) ")

        if new_m.log_Li < m.log_Li && !box_rflct
            process_report!(tuner, false)

            west_pos,_=box_move(m.pos,-rd,e.box)
            new_m=e.model_initλ(t, i, to_prior.(west_pos,e.priors), west_pos, -v′, e.obs, e.constants...)  #if east reflection fails, try west
            # println("newθ: $(to_prior.(west_pos,e.priors))")

            # println("west: $(new_m.log_Li-m.log_Li) ")

        end
        if new_m.log_Li < m.log_Li
            process_report!(tuner, false)
            south_pos,_=box_move(m.pos,-d,e.box)
            new_m=e.model_initλ(t, i, to_prior.(south_pos,e.priors),south_pos,-m.v, e.obs, e.constants...)  #if west reflection fails, try south
            # println("newθ: $(to_prior.(south_pos,e.priors))")

            # println("south: $(new_m.log_Li-m.log_Li) ")

        end
        new_m.log_Li < m.log_Li && (process_report!(tuner, false); m.v=-m.v) #if that fails reverse the particle velocity and wait for smaller τ
    end
    new_m.log_Li < m.log_Li ? (return m) : (process_report!(tuner, true); return new_m)
end

                function box_move(pos,d,box)
                    if !(all(box[:,1].<pos+d.<box[:,2]))
                        b=isapprox.(pos,box)
                        pos[b[:,1]].=box[b[:,1],1]
                        pos[b[:,2]].=box[b[:,2],2]
                        boundary_d=box.-pos
                        boundary_t=boundary_d./d
                        #boundary_t[b].=eps()
                        first_b_idx=findfirst(isequal(minimum(boundary_t[findall(x->x>=0.,boundary_t)])),boundary_t)
                        first_b_d=boundary_d[first_b_idx]
                        d=d.*(first_b_d/d[first_b_idx[1]])
                    end

                    return pos+d, d
                end

                function box_reflect(pos,box,v)
                    if !(all(box[:,1].<pos.<box[:,2]))
                        b=isapprox.(pos,box)
                        low_idxs=[v[i]<0&&b[i,1] for i in 1:length(v)]
                        hi_idxs=[v[i]>0&&b[i,2] for i in 1:length(v)]
                        v[low_idxs]=-v[low_idxs]
                        v[hi_idxs]=-v[hi_idxs]
                    end
                    return v
                end

                function boundary_norm(v, lhδ)
                    b=[lhδ/vi for vi in v]
                    b[b.==Inf].=prevfloat(Inf) #rectify invalid vals
                    b[b.==0.].=rand([nextfloat(0.),prevfloat(0.)])
                    b[b.==-Inf].=nextfloat(-Inf)  
                    return normalize(b)
                end

                function reflect(v, n, η)
                    v′=v-(2*n*dot(n,v))
                    η > 0. ? (return r_perturb(v′,η)) : (return v′)
                end

                function r_perturb(v′, η)
                    return v′*cos(η) + rand(MvNormal(length(v′),1.))*sin(η)
                end