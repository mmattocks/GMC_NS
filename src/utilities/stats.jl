function get_lognormal_params(desired_μ, desired_σ)
    μ²=desired_μ^2
    σ²=desired_σ^2
    ln_μ=log(μ²/sqrt(μ²+σ²))
    ln_σ=sqrt(log(1+(σ²/μ²)))
    return ln_μ, ln_σ
end

function marginals(d::NormalInverseGamma)
    m_T=MarginalTDist(d.mu,sqrt((d.shape*inv(d.v0))/d.scale),TDist(2*d.shape))
    m_Ga=InverseGamma(d.shape,d.scale)
    return m_T,m_Ga
end

function marginals(d::NormalGamma)
    m_T=MarginalTDist(d.mu,sqrt(d.rate/(d.shape*d.nu)),TDist(2*d.shape))
    m_Ga=Gamma(d.shape,1/d.rate)
    return m_T,m_Ga
end