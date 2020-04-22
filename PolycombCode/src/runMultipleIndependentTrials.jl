using Distributions # install by Pkg.add("Distributions") in julia
using Dates
using JLD
using HDF5
using DataFrames
using CSV
using Combinatorics
using Printf
using Distributed
using Random
using DelimitedFiles
using HypothesisTests

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
#### Run these 4 lines below when want to run code in julia terminal in Atom ####
configfile = "constants.jl"
indir = joinpath("..","input")
constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

include("utilities.jl")
include("types.jl")
include("individuals.jl")
include("population.jl")
include("measure.jl")
include("ipaNetworkFileConversion.jl")

#----------------------------------------------------------------
# Configure dataOutput folder for intact and broken polycomb data
#----------------------------------------------------------------
# Generate time stamp used to name folder to save results in
if RUNWITHANDWITHOUTPOLYCOMB
    timestamp = ARGS[1]
    outputVersion = if SAURABVERSION; "Saurabh"; else; "New"; end # If using Saurabh version or using new version where model environment explicitly
    outPcgStatus = if POLYCOMBFLAG; "Pcg"; else; "NoPcg"; end # If running model with Polycomb or without Polycomb
    dataIntactAndBrokenoutdir = ARGS[2]
    dataIntactAndBrokensimdir = joinpath(dataIntactAndBrokenoutdir, join([timestamp,outputVersion,outPcgStatus]))
    if !isdir(dataIntactAndBrokensimdir)
        mkdir(dataIntactAndBrokensimdir)
    end
else
    timestamp = gentimestamp()
    outputVersion = if SAURABVERSION; "Saurabh"; else; "New"; end # If using Saurabh version or using new version where model environment explicitly
    outPcgStatus = if POLYCOMBFLAG; "Pcg"; else; "NoPcg"; end # If running model with Polycomb or without Polycomb
    dataIntactAndBrokenoutdir = joinpath("..", "dataOutput")
    dataIntactAndBrokensimdir = joinpath(dataIntactAndBrokenoutdir, join([timestamp,outputVersion,outPcgStatus]))
    # save constants.jl file used for run of this code
    configpath = joinpath(dataIntactAndBrokensimdir, configfile)
    if !isdir(dataIntactAndBrokensimdir)
        mkdir(dataIntactAndBrokensimdir)
    end
    cp(constantsFile, configpath)
end
# Print out the folder data is being saved to:
print(dataIntactAndBrokensimdir)
#----------------------------------------------------------------

#----------------------------------------------------------------
# Initialize environmental and genetic robustness matrix to
# store results for each independent trial. As well as path
# length, number of different development state vectors, etc.
vecMeasurementsToTake = ["fitness","fitnessUnderEachEnv","develtypes","pcgStateTypes","pathlength","connectivity","robustness","envRobustness"]#ARGS[3]
data = Dict{String, DataFrame}()
for j = 1:DIFENVS
    for i=1:length(vecMeasurementsToTake)
        data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)] = DataFrame()
    end
end



#----------------------------------------------------------------
# Generate file path to the folder containing the .jld files for the
# independent populations previously generated using
# generateSavePopulations() function to run for each independent trial:
folderToGetIndependentPops = join(["savedPopsForIndTrialsWith", INDTRIALSFOLDER])
#----------------------------------------------------------------
# Start for loop for running each independent trial:
#----------------------------------------------------------------
for indTrial = 1:INDTRIALS
    if RANDSEEDFLAG == true
        Random.seed!(RANDSEEDNUM[indTrial]);
    end
    #----------------------------------------------------------------
    # Open initial populations composed of N individuals that were
    # previously saved in order to properly compare results between
    # different runs of runMultipleIndependentTrails.jl
    #----------------------------------------------------------------
    if SAVEDPOPSFLAG
        popFileName = join(["pop",indTrial,".jld"])
        popFolder = joinpath("../savedPopulationsForIndependentTrials",folderToGetIndependentPops)
        popFile = joinpath(popFolder,popFileName)
        pop1  = jldopen(popFile) do file
                    read(file, "Population")
                end
    else
        pop1 = genpop()
    end

    meas1 = genmeasure() # set the initial Measure type for each individual in the population
    # Dictory for all possible measurements:
    measureTypes = Dict("time"=>meas1.time, "fitness"=>meas1.fitness, "fitnessstd"=>meas1.fitnessstd, "fitnessUnderEachEnv"=>meas1.fitnessUnderEachEnv,
                    "fitnessUnderEachEnvStd"=>meas1.fitnessUnderEachEnvStd, "robustness"=>meas1.robustness, "robustnessstd"=>meas1.robustnessstd, "envRobustness"=>meas1.envRobustness,
                    "envRobustnessStd"=>meas1.envRobustnessStd, "pathlength"=>meas1.pathlength, "pathlengthstd"=>meas1.pathlengthstd,
                    "indtypes"=>meas1.indtypes, "inittypes"=>meas1.inittypes, "develtypes"=>meas1.develtypes, "pcgVecTypes"=>meas1.pcgVecTypes,
                    "pcgStateTypes"=>meas1.pcgStateTypes, "opttypes"=>meas1.opttypes, "connectivity"=>meas1.connectivity)
    #---------------------------------------------------------------------
    # Robustness Analysis
    #---------------------------------------------------------------------
    # First generate matrix with random different environment state vectors
    # that are only 4 environment changes away from initialized environment
    # state vector used during development and evolution, i.e. the founder's
    # environment state vector (for extended model).
    # OR generate matrix using initial state vector for Saurabh's version
    # Analysis for both intact and broken polycomb below use the same
    # different environment state vectors.
    # AND
    # Calculate genetic and environmental sensitivity of the founder to use to
    # measure genetic and environemtnal robustness for each individual
    if "robustness" in vecMeasurementsToTake
        founderGeneticRobust(pop1.founder)
        founderGeneticRobustness1 = pop1.founder.robustness # If evolving with multiple environments then this is vector containing founder robustness under each different environmental condition
    else
        founderGeneticRobustness1 = zeros(Float64,DIFENVS)
    end
    if "envRobustness" in vecMeasurementsToTake
        mutatedEnvMatrix1 = Array{Array{Int64,2},1}(undef,DIFENVS)
        for i = 1:DIFENVS
            mutatedEnvMatrix1[i] = generateMutatedEnvsMatrix(pop1.founder, i) # generate mutated
            # environemnt matrix based off of environment state vector (extended
            # version) OR initial state vector (saurabh version).
        end
        founderEnvRobust(pop1.founder)
        founderEnvRobustness1 = pop1.founder.envRobustness # will be vector of length DIFENVS  when evolving under multiple dif environmental conditions
    else
        founderEnvRobustness1 = zeros(Float64,DIFENVS)
        mutatedEnvMatrix1 = Array{Array{Int64,2},1}(undef,DIFENVS)
        for i = 1:DIFENVS
            mutatedEnvMatrix1[i] = Array{Float64,2}(undef,0,0)
        end
    end

    #------------------------------------------------------------------
    # Evolve population for GENS generations
    #------------------------------------------------------------------
    # Print out progress of running multiple independent trials:
    print("\nstart evolution for ", indTrial, " out of ", INDTRIALS, " independent trials at ", Dates.format(now(), "u dd, yy HH:MM:SS"))

    # Set the starting generation number to 1 before evolution of initial population begins
    measnum = 1 # initialize 1st measurement going to be taking every MEASPERIOD generations

    # avgIntactPolyResultsDuringEvolution = DataFrame(avgEnvRobustIntact = Float64[], avgPercentPliant = Float64[], avgPliantInact = Float64[], avgRobustInact = Float64[], avgStableEnvs = Float64[], avgNonStableEnvs = Float64[])
    for t=1:GENS
    # Running update function for total number of generations (GENS) --> evolution step

        update(pop1)
        if (mod(t-1, MEASPERIOD) == 0) | (t == GENS) # if the remainder of number of current generations/MEASPERIOD is zero OR t equals total number of generations (GENS)
            measure(pop1, meas1, t, measnum, mutatedEnvMatrix1, founderGeneticRobustness1, founderEnvRobustness1, vecMeasurementsToTake, measureTypes)
            measnum += 1
        end

    end
    for j = 1:DIFENVS
        for i=1:length(vecMeasurementsToTake)
            if length(getfield(meas1,Symbol(vecMeasurementsToTake[i]))[1,:]) > 1
                data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)][!, Symbol(string("trial",indTrial))] = getfield(meas1,Symbol(vecMeasurementsToTake[i]))[:,j]
            else
                data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)][!, Symbol(string("trial",indTrial))] = getfield(meas1,Symbol(vecMeasurementsToTake[i]))
            end
        end
    end
end

#----------------------------------------------------------------
# Save resulting meas (measurements taken every MEASPERIOD generations) for each independent trial
evolutionMeasurementsDir = joinpath(dataIntactAndBrokensimdir, "measurementsDuringEvolution")
if !isdir(evolutionMeasurementsDir)
    mkdir(evolutionMeasurementsDir)
end
# calculate and store the average and standard deviation environmental and genetic robustness for all the independent trials:
# could probably convert these for loop functions to using map or apply functions to help speed up code -ML 4/4/19
for j = 1:DIFENVS
    for i=1:length(vecMeasurementsToTake)
    # With DataFrames version that is older than v0.19.0 then usee colwise(), otherwise use [mean(col) for col in eachcol()]
    # To look at package version: pkgs = Pkg.installed(); pkgs["DataFrames"];
    # Use with DataFrames older than v0.19.0, so with Julia 0.7.0 or less
        # data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Matrix, data[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
        # data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Sem] = colwise(standardErrorOfMean, convert(DataFrame, transpose(convert(Matrix, data[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    # Use with DataFrames v0.19.0 or newer, so with Julia 1.0.0 or newer:
        data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)][!, :Average] = [mean(col) for col = eachcol(convert(DataFrame, transpose(convert(Matrix, data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)]))))]
        data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)][!, :Sem] = [standardErrorOfMean(col) for col = eachcol(convert(DataFrame, transpose(convert(Matrix, data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)]))))]
    end
end

for j = 1:DIFENVS
    averagesDuringEvolution = DataFrame()
    for i=1:length(vecMeasurementsToTake)
        averagesDuringEvolution[!, Symbol(string(vecMeasurementsToTake[i],"Avg","Env",j))] = data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)][!, :Average]
        averagesDuringEvolution[!, Symbol(string(vecMeasurementsToTake[i],"Sem","Env",j))] = data[string(vecMeasurementsToTake[i],"DuringEvolution","Env",j)][!, :Sem]
    end
    filetimestampAndVersions = join([timestamp,outputVersion,outPcgStatus])
    filenameForDifEnvs = string("avgsDuringEvolutionMeasurements",j,"Env.csv")
    CSV.write(joinpath(evolutionMeasurementsDir, filenameForDifEnvs), averagesDuringEvolution)
end

#-------------------------------------------------------------------------
# print out the name of the folder saved in DataOutput folder for this
# particular run of multiple independent trials.
print("\nData saved to:")
print("\n",dataIntactAndBrokensimdir)
print("\n===========================================\n")
