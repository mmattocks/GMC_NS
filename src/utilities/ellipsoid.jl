# Represents an N-dimensional ellipsoid
struct Ellipsoid
    ctr::Vector{Float64}  # center coordinates
    cov::Array{Float64, 2}
    icov::Array{Float64, 2}  # inverse of cov
    vol::Float64
end

# Draw a random point from within a unit N-ball
function randnball(ndim)
    z = randn(ndim)
    r2 = 0.
    for i=1:ndim
        r2 += z[i]*z[i]
    end
    factor = rand()^(1. /ndim) / sqrt(r2)
    for i=1:ndim
        z[i] *= factor
    end
    return z
end

# proportionality constant depending on dimension
# for n even:      (2pi)^(n    /2) / (2 * 4 * ... * n)
# for n odd :  2 * (2pi)^((n-1)/2) / (1 * 3 * ... * n)
function nball_vol_factor(ndim::Int)
    if ndim % 2 == 0
        c = 1.
        for i=2:2:ndim
            c *= 2pi / i
        end
        return c
    else
        c = 2.
        for i = 3:2:ndim
            c *= 2pi / i
        end
        return c
    end
end

function ellipsoid_volume(scaled_cov::Matrix{Float64})
    ndim = size(scaled_cov, 1)
    return nball_vol_factor(ndim) * sqrt(det(scaled_cov))
end

# find the bounding ellipsoid of points x where 
function bounding_ellipsoid(x::Matrix{Float64}, enlarge=1.0)

    ndim, npoints = size(x)

    ctr = mean(x, dims=2)[:, 1]
    delta = x .- ctr
    cov = unscaled_covzm(delta, 2)
    icov = inv(cov)

    # Calculate expansion factor necessary to bound each point.
    # This finds the maximum of (delta_i' * icov * delta_i)
    fmax = -Inf
    for k in 1:npoints
        f = 0.0
        for j=1:ndim
            for i=1:ndim
                f += icov[i, j] * delta[i, k] * delta[j, k]
            end
        end
        fmax = max(fmax, f)
    end

    fmax *= enlarge
    cov .*= fmax
    icov .*= 1. /fmax
    vol = ellipsoid_volume(cov)

    return Ellipsoid(ctr, cov, icov, vol)
end

function sample_ellipsoid(ell::Ellipsoid)
    ndim = length(ell.ctr)

    # Get scaled eigenvectors (in columns): vs[:,i] is the i-th eigenvector.
    f = eigen(ell.cov)
    v, w = f.vectors, f.values
    for j=1:ndim
        tmp = sqrt(abs(w[j]))
        for i=1:ndim
            v[i, j] *= tmp
        end
    end

    return v*randnball(ndim) + ell.ctr
end