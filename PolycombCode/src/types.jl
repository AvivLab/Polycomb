abstract type Individuals end
#

mutable struct Individual{T} <: Individuals
    network::Matrix{T} # W matrix with GenesXGenes representing genes interactions
    initstate::Vector{Float64} # Sw1 vector: initial state vector
    EnvIndInter::Matrix{T} # W x E with GenesxEnivornments representing environment gene interactions
    EnvState::Matrix{T} # Se1 vector; environment state vectors where each column is different environment state vector, like 2nd column would be Se2 vector; environment state vector/input that is 70% different than EnvState
    develstate::Matrix{T} # output of individual state vector of Genes length, after maxConv iterations with entries that are expression level of each gene
    optstate::Matrix{T} # optimum state - stable founder's development state vector (found using generate founder function at beginning of simulation)
    stable::BitArray{} # True if Individual is stable determined if current initial state is equal to the updated initial state vector (so last iteration equal to previous iteration state vector)
    fitness::Float64 # total fitness of individual under all different environmental states
    fitnessUnderEachEnv::Vector{Float64} # fitness's of an individual under each different environmental condition separately
    pathlength::Vector{Int64} # number iterations needed to reach stability in E1
    polycombvec::Matrix{T} # theta vector if only one PRC of polycomb susceptible genes: 0 if gene not susceptible to polycomb and 1 if gene is susceptible to polycomb; matrix when have multiple PRCs where column is each PRC and each row is that PRC's effect on a given gene
    polycombstate::Matrix{T} # after TIMECRIT what polycomb susceptible genes are expressed or supressed: 0 if susceptible gene is suppressed (meaning polycomb response element active) in given env, and 1 if susceptible gene is expressed (meaning polycomb response element is inactive) in given env. So matrix where rows are genes are columns and envs
    robustness::Vector{Float64} # measure of genetic robustness to single mutations
    envRobustness::Vector{Float64} # measure of robustness to changing environments
end


mutable struct Population{I <: Individuals}
    individuals::Vector{I}
    founder::I
end


mutable struct Measure
    time::Vector{Int64} # generation take measurement(s)
    fitness::Vector{Float64} # total fitness of individual under multiple different environmental conditions
    fitnessstd::Vector{Float64}
    fitnessUnderEachEnv::Matrix{Float64} # fitness of individual under each environment condition separately
    fitnessUnderEachEnvStd::Matrix{Float64}
    robustness::Matrix{Float64}
    robustnessstd::Matrix{Float64}
    envRobustness::Matrix{Float64}
    envRobustnessStd::Matrix{Float64}
    pathlength::Matrix{Float64}
    pathlengthstd::Matrix{Float64}
    indtypes::Vector{Float64} # number of different "networks" (W matrices)
    inittypes::Vector{Float64} # number of different initialStateVecs
    develtypes::Matrix{Float64} # number of different develstate vectors
    pcgVecTypes::Vector{Float64}
    pcgStateTypes::Matrix{Float64}
    opttypes::Matrix{Float64}
    connectivity::Vector{Float64}
end

#const Time = Int64
