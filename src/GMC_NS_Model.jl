"""
    GMC_NS_Model_Record

Abstract supertype for records of serialised 'GMC_NS_Model's. Subtyped record structs must implement the following fields as their interface with 'GMC_NS':

    struct Minimal_Record
        trajectory::Integer #id assigned at construction, identifies trajectory
        i::Integer #id integer assigned at construction, identifies position on trajectory
        pos::Vector{<:AbstractFloat} #hypersphere position of model-particle
        path::String #model is serialized to ensembledir/trajectory.i
        log_Li::AbstractFloat #log likelihood calculated by 'GMC_NS_Model' constructor
    end
"""
abstract type GMC_NS_Model_Record end

"""
    GMC_NS_Model_Record

Abstract supertype for models in 'GMC_NS_Ensemble's. Subtyped model structs must implement the following fields as their interface with 'GMC_NS':

    mutable struct Minimal_Model
        trajectory::Integer #id assigned at construction, identifies trajectory
        i::Integer #id integer assigned at construction, identifies position on trajectory
        
        θ::Vector{<:AbstractFloat} #model's parameter Vector
        log_Li<:AbstractFloat #model's likelihood, calculated by constructor

        pos::Vector{<:AbstractFloat}
        v::Vector{<:AbstractFloat} #model's velocity in hypersphere coords
    end
"""
abstract type GMC_NS_Model end