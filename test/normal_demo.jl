using GMC_NS, Distributions, ConjugatePriors

sample_dist=Normal(100.,5.)
samples=rand(sample_dist, 10000)

prior=NormalGamma(300.,2e-3,1.,10.) #a lousy prior

box=[0. 1000.
     1/1000. 1/eps()]

gmc=GMC_DEFAULTS
gmc[1]=5
gmc[2]=eps()

e=Normal_Ensemble("NE_test", 10, samples, prior, box, gmc...)

uds=Vector{Vector{Function}}([[convergence_display],[evidence_display],[info_display]])

lds=Vector{Vector{Function}}([[model_obs_display]])

converge_ensemble!(e,backup=(true,100),upper_displays=uds, lower_displays=lds, disp_rot_its=1000)