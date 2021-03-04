## Galilean Monte Carlo Nested Sampling
[![Build Status](https://travis-ci.org/mmattocks/GMC_NS.jl.svg?branch=master)](https://travis-ci.org/mmattocks/GMC_NS.jl)
[![codecov](https://codecov.io/gh/mmattocks/GMC_NS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmattocks/GMC_NS.jl)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

Implements Galilean Monte Carlo Nested Sampling in pure Julia, intended for arbitrary models with ordinary continuous parameter spaces.

Skilling, John. “Galilean and Hamiltonian Monte Carlo.” In Proceedings, 33:8. Garching, Germany: MDPI, 2019.

##Installation
To install, add the github address in the Julia package manager (shortcut "]").
'''
julia> ]add https://github.com/mmattocks
'''

##Example: Solving the Eggbox Problem

The eggbox problem uses the likelihood function 'log_lh=(2 + cos(5π * x1) * cos(5π * x2))^5', where x1 and x2 are the two parameters of the model, in the range 0:1. This gives a surface with 16 evenly spaced minima, the "eggbox". To solve it with GMC_NS.jl:

'''
using GMC_NS, Distributions

#define the priors
prior=Uniform(0,1)
priors=[prior,prior]

#define the sampling box
box=[0. 1.
     0. 1.]

#use the package's default settings
gmc=GMC_DEFAULTS
gmc[1]=160 #change the minimum number of live particles to 160
gmc[2]=eps() #change the timestep at which particles are killed to be very small

#instantiate an ensemble with 1000 model-particles, given the priors, box, and GMC settings
e=Eggbox_Ensemble("Eggbox_test", 1000, priors, box, gmc...)

#define some displays to monitor the algorithm on-line
uds=Vector{Vector{Function}}([[evidence_display],[liwi_display],[convergence_display]])
lds=Vector{Vector{Function}}([[ensemble_display]])

#perform GMC-NS steps until the ensemble is converged
converge_ensemble!(e,backup=(true,100),upper_displays=uds,lower_displays=lds, disp_rot_its=1000, converge_factor=1e-6)
'''

Output:
![Eggbox output](https://github.com/mmattocks/GMC_NS.jl/gmcns.gif)