using GMC_NS, Distributions, ConjugatePriors, ProgressMeter, Serialization

import GMC_NS:Normal_Model,Normal_Ensemble

sample_dist=Normal(100.,5.)
samples=rand(sample_dist, 10000)

prior=NormalGamma(300,1,1,.1) #a lousy prior

e=Normal_Ensemble("NE_test", 300, samples, [prior], GMC_DEFAULTS...)

uds=Vector{Vector{Function}}([[tuning_display],[convergence_display],[evidence_display],[info_display]])

lds=Vector{Vector{Function}}([[model_obs_display]])

converge_ensemble!(e,backup=(true,100),upper_displays=uds, lower_displays=lds, disp_rot_its=1000)