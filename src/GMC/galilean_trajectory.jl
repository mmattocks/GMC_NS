"""
    galilean_trajectory_sample(e, m, τ)

Given 'm <: GMC_NS_Model', 'e <: GMC_NS_Ensemble', and time step τ, attempt to find new model obeying 'new_m.log_Li > e.contour', first by permuting 'm.θ' by the distance given by velocity 'm.v' and time step 'τ', with 'm' retaining particle velocity v:

''[θ,v] -> [θ′=θ+τv, v]''

 If this fails, the rejected particle "reflects" off the contour boundary as represented by the gradient normal 'n':

 ''[θ,v] -> [θ′′=θ+τv+τv′, v′=v-2n(n∘v)]''

If 'e.GMC_PERTURB' is a positive value, this value (η below) is used with a newly sampled isotropic vector to calculate a perturbed reflection vector 'v′′':

''v′′=v′cosη+rsinη, r=Normal(0,I)''
 
Finally, if reflection fails to produce a new particle inside the contour, invert the particle's velocity and return it to the ensemble, preserving detailed balance:

''[θ,v] -> [θ, v′=-v]''

"""
function galilean_trajectory_sample!(m, e, τ)
    d=τ*m.v #distance vector for given timestep
    new_m=e.model_initλ(0, m.id, m.θ+d, m.v, e.constants...) #try to proceed along distance vector
    if new_m.log_Li <= e.contour #if that didnt work
        lhδ=m.log_Li-new_m.log_Li #get the lhδ between the two points
        n=boundary_norm(d,lhδ) #get the gradient normal for the distance and lhδ
        v′=reflect(m.v,n,e.GMC_reflect_η) #get the reflected velocity vector off the boundary
        rd=τ*v′ #reflected distance given the timestep
        new_m=e.model_initλ(0, m.id, m.θ+d+rd, v′, e.constants...)

        if new_m.log_Li <= e.contour #if reflection fails
            m.v=-m.v #reverse the particle's velocity
        end
    end
    new_m.log_lh <= e.contour ? (return m) : (return new_m)
end

                function boundary_norm(v, lhδ)
                    return normalize([lhδ/vi for vi in v])
                end

                function reflect(v, n, η)
                    v′=v-(2*n*dot(n,v))
                    η > 0. ? (return r_perturb(v′,η)) : (return v′)
                end

                function r_perturb(v′, η)
                    return v′*cos(η) + rand(MvNormal(length(v′),1.))*sin(η)
                end
