"""
    GMC_NS_Model_Record

Abstract supertype for records of serialised 'GMC_NS_Model's. Subtyped record structs must implement the following fields as their interface with 'GMC_NS':

    struct Minimal_Record
        path::String #model is serialized here
        log_Li::Float64 #log likelihood calculated by 'GMC_NS_Model' constructor
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
        v::Vector{<:AbstractFloat} #model's velocity in parameter space
        log_Li<:AbstractFloat #model's likelihood, calculated by constructor
    end
"""
abstract type GMC_NS_Model end