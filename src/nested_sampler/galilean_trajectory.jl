function galilean_trajectory_sample(e, m, τ, contour)
    new_θ=m.θ+τ*m.v
    new_m=e.model_initλ(new_θ, e.constants...)
    if new_m.log_lh <= contour
        lhδ=m.log_lh-new_m.log_lh
        new_θ+=reflect(τ*m.v,boundary_norm(v,lhδ))
        new_m=e.model_initλ(new_θ, e.constants...)

        if new_m.log_lh <= contour
            new_θ=m.θ+τ*(-m.v)
            new_m=e.model_initλ(new_θ, e.constants...)
        end
    end

    new_m.log_lh <= contour ? (return nothing) : (return new_m)
end

function boundary_norm(v, lhδ)
    return normalize([lhδ/vi for vi in v])
end

function reflect(v, n; perturb=true)
    peturb ? (return r_perturb(v-(2*n*dot(n,v)))) : (return v-(2*n*dot(n,v)))
end

function r_perturb(v′; θ=.001)
    return v′*cos(θ) + rand(MvNormal(length(v′),1.))*sin(θ)
end
