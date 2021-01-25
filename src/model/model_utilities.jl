function construct_reflection(m, i, r_v)
    params=[getfield(m,pn) for pn in propertynames(m)]
    params[2]=i
    params[6]=r_v
    return typeof(m)(params...)
end
