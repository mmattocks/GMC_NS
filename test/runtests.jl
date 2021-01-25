@info "Loading test packages..."
using GMC_NS, Test, Distributions, Serialization
import GMC_NS:box_move, box_reflect, boundary_norm, reflect, r_perturb
import LinearAlgebra:normalize

@info "Beginning tests..."
using Random
Random.seed!()

include("galilean_unit_tests.jl")
include("ensemble_tests.jl")

@info "Tests complete!"
