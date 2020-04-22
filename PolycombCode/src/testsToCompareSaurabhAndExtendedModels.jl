using Base.Test
using Distributions
using Base.Dates
using DataFrames
#using GR
using Plots
using Combinatorics


configfile = "constants.jl"
indir = joinpath("..","input")
constantsFile = joinpath(indir,configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

include("types.jl")
include("individuals.jl")
include("population.jl")
include("measure.jl")

################################################################################
# Test environmental robustness behavior of population after
# evolution --> Saurabh version compared to extended version
################################################################################
## Test environmentalRobustness() function
SAURABVERSION = true
const WMAT = Matrix{Float64}(0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
saurabhPop = genpop()
mutatedEnvMatrix = generateMutatedEnvsMatrix(saurabhPop.founder)
founderEnvRobustness = founderEnvRobust(saurabhPop.founder, mutatedEnvMatrix)
individualRobustnessBeforeEvol = pmap(i -> founderEnvRobust(i, mutatedEnvMatrix), saurabhPop.individuals)
@test Float32(mean(individualRobustnessBeforeEvol)) == Float32(founderEnvRobustness)
for i = 1:20
    update(saurabhPop)
end
sensitivitySaurabh = map(i -> founderEnvRobust(saurabhPop.individuals[i],mutatedEnvMatrix), 1:N)
robustnessCalculatedSaurabh = founderEnvRobustness - mean(sensitivitySaurabh)
robustnessInCodeSaurabh = map(i -> environmentalRobustness(saurabhPop.individuals[i],mutatedEnvMatrix,founderEnvRobustness), 1:N)
abs(mean(robustnessInCodeSaurabh) - robustnessCalculatedSaurabh)


## Compare Saurabh's C++ code version of calculating robustness: founder.develstate - perturbed.develstate instead of me.develstate - perturbed.develstate
SAURABVERSION = true
const WMAT = Matrix{Float64}(0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
saurabhPop = genpop()
mutatedEnvMatrix = generateMutatedEnvsMatrix(saurabhPop.founder)
founderEnvRobustness = founderEnvRobust(saurabhPop.founder, mutatedEnvMatrix) # same for us and Saurabh's C++ code
for i = 1:20
    update(saurabhPop)
end
differencesInMethods = map(i -> environmentalRobustnessSaurabhCcodeVersion(saurabhPop.individuals[i],mutatedEnvMatrix,founderEnvRobustness,saurabhPop.founder) - environmentalRobustness(saurabhPop.individuals[i],mutatedEnvMatrix,founderEnvRobustness), 1:N)

# Looks like typically Saurabh's C++ version causes env robustness to be less
# than our version...This makes sense because as evolution goes on each
# individual's develstate can get farther and farther from the founders
# develstate because not strong selection (SELSTR = 1) (fitness decreases
# during evolution), so sensitivity will increase with evolution and thus
# robustness will decrease...This seems kinda like cheating to get better
# results... -ML 4/4/19
t =1
## Function to run evolution to find environment robustness for plotting behavior:
function runEvolutionToFindEnvRobustness(pop::Population,sensivityflag::Bool)
    vecMeasurementsToTake = ["envRobustness"]
    meas = genmeasure() # set the initial Measure type for each individual in the population
    measureTypes = Dict("time"=>meas.time, "fitness"=>meas.fitness, "fitnessstd"=>meas.fitnessstd,
                    "robustness"=>meas.robustness, "robustnessstd"=>meas.robustnessstd, "envRobustness"=>meas.envRobustness,
                    "envRobustnessStd"=>meas.envRobustnessStd, "pathlength"=>meas.pathlength, "pathlengthstd"=>meas.pathlengthstd,
                    "indtypes"=>meas.indtypes, "inittypes"=>meas.inittypes, "develtypes"=>meas.develtypes, "pcgVecTypes"=>meas.pcgVecTypes,
                    "pcgStateTypes"=>meas.pcgStateTypes, "opttypes"=>meas.opttypes, "phenotypicDist"=>meas.phenotypicDist,
                    "phenotypicDistStd"=>meas.phenotypicDistStd, "connectivity"=>meas.connectivity)
    # Set the starting generation number to 1 before evolution of initial population begins
    measnum = 1 # initialize 1st measurement going to be taking every MEASPERIOD generations
    # Calculate genetic and environmental sensitivity of the founder to use to
    # measure genetic and environemtnal robustness for each individual
    founderGeneticRobustness = 0.
    mutatedEnvMatrix = generateMutatedEnvsMatrix(pop.founder)
    founderEnvRobustness = founderEnvRobust(pop.founder, mutatedEnvMatrix)
    # avgIntactPolyResultsDuringEvolution = DataFrame(avgEnvRobustIntact = Float64[], avgPercentPliant = Float64[], avgPliantInact = Float64[], avgRobustInact = Float64[], avgStableEnvs = Float64[], avgNonStableEnvs = Float64[])
    for t=1:GENS
    # Running update function for total number of generations (GENS) --> evolution step
        update(pop)
        if (mod(t-1, MEASPERIOD) == 0) | (t == GENS) # if the remainder of number of current generations/MEASPERIOD is zero OR t equals total number of generations (GENS)
            if sensivityflag
                envsensvect = pmap(i -> environmentalSensitivity(i, mutatedEnvMatrix), pop.individuals)
                meas.envRobustness[measnum] = mean(envsensvect)
            else
                envrobustvect = pmap(i -> environmentalRobustness(i, mutatedEnvMatrix, founderEnvRobustness), pop.individuals)
                meas.envRobustness[measnum] = mean(envrobustvect)
            end
            measnum += 1
        end
    end
    return meas.envRobustness
end
## Investigate using Saurabh's C++ code version when run evolution to find env Robustness:
function runEvolutionToGetMeasurements(pop::Population,vecMeasurementsToTake::Vector)
    meas = genmeasure()
    measCcode = genmeasure() # set the initial Measure type for each individual in the population
    measureTypes = Dict("time"=>meas.time, "fitness"=>meas.fitness, "fitnessstd"=>meas.fitnessstd,
                    "robustness"=>meas.robustness, "robustnessstd"=>meas.robustnessstd, "envRobustness"=>meas.envRobustness,
                    "envRobustnessStd"=>meas.envRobustnessStd, "pathlength"=>meas.pathlength, "pathlengthstd"=>meas.pathlengthstd,
                    "indtypes"=>meas.indtypes, "inittypes"=>meas.inittypes, "develtypes"=>meas.develtypes, "pcgVecTypes"=>meas.pcgVecTypes,
                    "pcgStateTypes"=>meas.pcgStateTypes, "opttypes"=>meas.opttypes, "connectivity"=>meas.connectivity)

    measureTypesCcode = Dict("time"=>measCcode.time, "fitness"=>measCcode.fitness, "fitnessstd"=>measCcode.fitnessstd,
                    "robustness"=>measCcode.robustness, "robustnessstd"=>measCcode.robustnessstd, "envRobustness"=>measCcode.envRobustness,
                    "envRobustnessStd"=>measCcode.envRobustnessStd, "pathlength"=>measCcode.pathlength, "pathlengthstd"=>measCcode.pathlengthstd,
                    "indtypes"=>measCcode.indtypes, "inittypes"=>measCcode.inittypes, "develtypes"=>measCcode.develtypes, "pcgVecTypes"=>measCcode.pcgVecTypes,
                    "pcgStateTypes"=>measCcode.pcgStateTypes, "opttypes"=>measCcode.opttypes, "connectivity"=>measCcode.connectivity)
    # Set the starting generation number to 1 before evolution of initial population begins
    measnum = 1 # initialize 1st measurement going to be taking every MEASPERIOD generations
    if "robustness" in vecMeasurementsToTake
        founderGeneticRobustness = founderGeneticRobust(pop.founder)
    else
        founderGeneticRobustness = 0.
    end
    if "envRobustness" in vecMeasurementsToTake
        mutatedEnvMatrix = generateMutatedEnvsMatrix(pop.founder) # generate mutated
        # environemnt matrix based off of environment state vector (extended
        # version) OR initial state vector (saurabh version).
        founderEnvRobustness = founderEnvRobust(pop.founder, mutatedEnvMatrix)
    else
        founderEnvRobustness = 0.
        mutatedEnvMatrix = Array{Float64,2}[]
    end
    for t=1:GENS
    # Running update function for total number of generations (GENS) --> evolution step
        update(pop)
        if (mod(t-1, MEASPERIOD) == 0) | (t == GENS) # if the remainder of number of current generations/MEASPERIOD is zero OR t equals total number of generations (GENS)
            measure(pop, meas, t, measnum, mutatedEnvMatrix, founderGeneticRobustness, founderEnvRobustness, vecMeasurementsToTake, measureTypes)
            measureUsingSaurabhCcodeVersionForEnvRobust(pop, measCcode, t, measnum, mutatedEnvMatrix, founderGeneticRobustness, founderEnvRobustness, vecMeasurementsToTake, measureTypesCcode)
            measnum += 1
        end
    end
    return meas, measCcode
end
vecMeasurementsToTake = ["envRobustness","robustness","fitness","develtypes","pcgStateTypes","pathlength"]
data = Dict{String, DataFrame}()
dataCcode = Dict{String, DataFrame}()
dataPcG = Dict{String, DataFrame}()
dataCcodePcG = Dict{String, DataFrame}()
dataStrongSel = Dict{String, DataFrame}()
dataCcodeStrongSel = Dict{String, DataFrame}()
dataStrongSelPcG = Dict{String, DataFrame}()
dataCcodeStrongSelPcG = Dict{String, DataFrame}()
for i=1:length(vecMeasurementsToTake)
    data[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataCcode[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataPcG[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataCcodePcG[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataCcodeStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataCcodeStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
end
for indTrial = 1:100
    print("indTrial = ", indTrial)
    SELSTR = 1.0
    POLYCOMBFLAG = false
    TIMECRIT = 0
    include("individuals.jl")
    saurabhPop = genpop()
    saurabhPopPcG = deepcopy(saurabhPop)
    saurabhPopStrongSel = deepcopy(saurabhPop)
    saurabhPopStrongSelPcG = deepcopy(saurabhPop)
    meas, measCcode = runEvolutionToGetMeasurements(saurabhPop, vecMeasurementsToTake)

    SELSTR = 0.2
    include("individuals.jl")
    measStrongSel, measCcodeStrongSel = runEvolutionToGetMeasurements(saurabhPopStrongSel, vecMeasurementsToTake)

    POLYCOMBFLAG = true
    TIMECRIT = 1
    include("individuals.jl")
    measStrongSelPcG, measCcodeStrongSelPcG = runEvolutionToGetMeasurements(saurabhPopStrongSelPcG, vecMeasurementsToTake)

    SELSTR = 1.0
    include("individuals.jl")
    measPcG, measCcodePcG = runEvolutionToGetMeasurements(saurabhPopPcG, vecMeasurementsToTake)

    for i=1:length(vecMeasurementsToTake)
        data[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(meas,Symbol(vecMeasurementsToTake[i]))
        dataCcode[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measCcode,Symbol(vecMeasurementsToTake[i]))
        dataPcG[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measPcG,Symbol(vecMeasurementsToTake[i]))
        dataCcodePcG[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measCcodePcG,Symbol(vecMeasurementsToTake[i]))
        dataStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measStrongSelPcG,Symbol(vecMeasurementsToTake[i]))
        dataCcodeStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measCcodeStrongSelPcG,Symbol(vecMeasurementsToTake[i]))
        dataStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measStrongSel,Symbol(vecMeasurementsToTake[i]))
        dataCcodeStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(measCcodeStrongSel,Symbol(vecMeasurementsToTake[i]))
    end
end
# calculate and store the average and standard deviation environmental and genetic robustness for all the independent trials:
# could probably convert these for loop functions to using map or apply functions to help speed up code -ML 4/4/19
for i=1:length(vecMeasurementsToTake)
    data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, data[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataCcode[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataCcode[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataPcG[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataPcG[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataCcodePcG[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataCcodePcG[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataCcodeStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataCcodeStrongSelPcG[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    dataCcodeStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Array, dataCcodeStrongSel[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
end
# Plot average hamming distance of develstates compared to founder develstate (aka optstate):
# "envRobustness","robustness","fitness","develtypes","pcgStateTypes","pathlength"
measurementName = "pathlength"
measurementToPlotString = string(measurementName,"DuringEvolution")
if measurementName == "envRobustness"
    saurabhPlot = plot(Vector(1:21),Array(data[measurementToPlotString][:Average]),title="No PcG, SELSTR = 1",legend=false)
    saurabhPlotStrongSel = plot(Vector(1:21),Array(dataStrongSel[measurementToPlotString][:Average]),title="No PcG, SELSTR = 0.2",legend=false)
    saurabhPlotStrongSelPcG = plot(Vector(1:21),Array(dataStrongSelPcG[measurementToPlotString][:Average]),title="PcG, SELSTR = 0.2",legend=false)
    saurabhPlotPcG = plot(Vector(1:21),Array(dataPcG[measurementToPlotString][:Average]),title="PcG, SELSTR = 1",legend=false)
    saurabhPlotCcode = plot(Vector(1:21),Array(dataCcode[measurementToPlotString][:Average]),title="C++, No PcG, SELSTR = 1",legend=false)
    saurabhPlotCcodeStrongSel = plot(Vector(1:21),Array(dataCcodeStrongSel[measurementToPlotString][:Average]),title="C++, No PcG, SELSTR = 0.2",legend=false)
    saurabhPlotCcodeStrongSelPcG = plot(Vector(1:21),Array(dataCcodeStrongSelPcG[measurementToPlotString][:Average]),title="C++, PcG, SELSTR = 0.2",legend=false)
    saurabhPlotCcodePcG = plot(Vector(1:21),Array(dataCcodePcG[measurementToPlotString][:Average]),title="C++, PcG, SELSTR = 1",legend=false)
    annotate!([(20,0.0, text(measurementToPlotString, 16, :red, :center))])
    plot(saurabhPlot,saurabhPlotCcode,saurabhPlotPcG,saurabhPlotCcodePcG,saurabhPlotStrongSel,saurabhPlotCcodeStrongSel,
        saurabhPlotStrongSelPcG,saurabhPlotCcodeStrongSelPcG,layout=(4,2))
else
    saurabhPlot = plot(Vector(1:21),Array(data[measurementToPlotString][:Average]),title="No PcG, SELSTR = 1",legend=false)
    saurabhPlotStrongSel = plot(Vector(1:21),Array(dataStrongSel[measurementToPlotString][:Average]),title="No PcG, SELSTR = 0.2",legend=false)
    saurabhPlotStrongSelPcG = plot(Vector(1:21),Array(dataStrongSelPcG[measurementToPlotString][:Average]),title="PcG, SELSTR = 0.2",legend=false)
    saurabhPlotPcG = plot(Vector(1:21),Array(dataPcG[measurementToPlotString][:Average]),title="PcG, SELSTR = 1",legend=false)
    annotate!([(20,0.0, text(measurementToPlotString, 16, :red, :center))])
    plot(saurabhPlot,saurabhPlotPcG,saurabhPlotStrongSel,saurabhPlotStrongSelPcG,layout=(2,2))
end


# Plot averagees to compare selection pressures:
saurabhStrongAvgPlot = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhStrongSel[:envRobustAvg]),title="Saurabh strong sel Avg Robust",legend=false)
saurabhAvgPlot = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabh[:envRobustAvg]),title="Saurabh Avg Robust",legend=false)
saurabhStrongAvgPlotCcode = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhCcodeStrongSel),title="Saurabh Ccode strong sel Avg Robust",legend=false)
saurabhAvgPlotCcode = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhCcode),title="Saurabh Ccode Avg Robust",legend=false)
plot(saurabhStrongAvgPlot,saurabhStrongAvgPlotCcode,saurabhAvgPlot,saurabhAvgPlotCcode,layout=(2,2))
# Plot extended:
extendedPlotCcode = plot(Vector(1:21),Array(envRobustnessDuringEvolutionExtendedSaurabhCcode),title="Extended Model C code Env Robustness")
saurabhPlotCcode = plot(Vector(1:21),Array(envRobustnessDuringEvolutionSaurabhCcode),title="Saurabh Model C code Env Robustness")
plot(extendedPlotCcode,saurabhPlotCcode,layout=(2,1))
# Plot our code results:
extendedPlot = plot(Vector(1:21),Array(envRobustnessDuringEvolutionExtended),ylims=(-0.05,0.05),title="Extended Model Env Robustness")
saurabhPlot = plot(Vector(1:21),Array(envRobustnessDuringEvolutionSaurabh),title="Saurabh Model Env Robustness")
plot(extendedPlot,saurabhPlot,layout=(2,1))
# Plot averages for env robustness:
extendedAvgPlot = plot(Vector(1:21),Array(averagesDuringEvolutionExtended),title="Extended Avg Robust")
saurabhAvgPlot = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabh[:envRobustAvg]),title="Saurabh Avg Robust")
extendedAvgPlotCcode = plot(Vector(1:21),Array(averagesDuringEvolutionExtendedSaurabhCcode),title="Extended Ccode Avg Robust")
saurabhAvgPlotCcode = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhCcode),title="Saurabh Ccode Avg Robust")
plot(extendedAvgPlotCcode,saurabhAvgPlotCcode,extendedAvgPlot,saurabhAvgPlot,layout=(2,2))
# Plot averages for env robustness without polycomb:
extendedAvgPlotWoPcg = plot(Vector(1:21),Array(averagesDuringEvolutionExtendedWoPcg),title="Extended Avg WoPcg Robust")
saurabhAvgPlotWoPcg = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhWoPcg),title="Saurabh Avg WoPcg Robust")
extendedAvgPlotCcodeWoPcg = plot(Vector(1:21),Array(averagesDuringEvolutionExtendedSaurabhCcodeWoPcg),title="Extended Ccode Avg WoPcg Robust")
saurabhAvgPlotCcodeWoPcg = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhCcodeWoPcg),title="Saurabh Ccode Avg WoPcg Robust")
plot(extendedAvgPlotCcodeWoPcg,saurabhAvgPlotCcodeWoPcg,extendedAvgPlotWoPcg,saurabhAvgPlotWoPcg,layout=(2,2))
# Plot averages for environmental sensitivity
# extendedAvgPlotSens = plot(Vector(1:21),Array(averagesDuringEvolutionExtendedSensitiv),title="Extended Avg Sens")
# saurabhAvgPlotSens = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhSensitiv),title="Saurabh Avg Sens")
# extendedAvgPlotCcodeSens = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhCcodeSensitiv),title="Extended Ccode Avg Sens")
# saurabhAvgPlotCcodeSens = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabhCcodeSensitiv),title="Saurabh Ccode Avg Sens")
# plot(extendedAvgPlotSens,saurabhAvgPlotSens,extendedAvgPlotCcodeSens,saurabhAvgPlotCcodeSens,layout=(2,2))
##

# envRobustnessDuringEvolutionSaurabh = DataFrame()
# envRobustnessDuringEvolutionExtended = DataFrame()
# envRobustnessDuringEvolutionExtendedButTreatedLikeSaurabh = DataFrame()
# for indTrial = 1:2
#     SAURABVERSION = true
#     const WMAT = Matrix{Float64}(0,0)
#     const INP = (Float64)[]
#     const OPT = (Float64)[]
#     include("types.jl")
#     include("individuals.jl")
#     include("population.jl")
#     saurabhPop = genpop()
#     envRobustnessDuringEvolutionSaurabh[Symbol(string("trial",indTrial))] = runEvolutionToFindEnvRobustness(saurabhPop)
#
#     SAURABVERSION = false
#     const WMAT = saurabhPop.founder.network
#     const INP = saurabhPop.founder.initstate1
#     const OPT = (Float64)[]
#     include("types.jl")
#     include("individuals.jl")
#     include("population.jl")
#     extendedPop = genpop()
#     envRobustnessDuringEvolutionExtended[Symbol(string("trial",indTrial))] = runEvolutionToFindEnvRobustness(extendedPop)
#
#     SAURABVERSION = true
#     include("types.jl")
#     include("individuals.jl")
#     include("population.jl")
#     envRobustnessDuringEvolutionExtendedButTreatedLikeSaurabh[Symbol(string("trial",indTrial))] = runEvolutionToFindEnvRobustness(extendedPop)
# end
# averagesDuringEvolutionExtended = DataFrame(envRobustAvg = colwise(mean, convert(DataFrame, transpose(convert(Array, envRobustnessDuringEvolutionExtended)))))
# averagesDuringEvolutionSaurabh = DataFrame(envRobustAvg = colwise(mean, convert(DataFrame, transpose(convert(Array, envRobustnessDuringEvolutionSaurabh)))))
# averagesDuringEvolutionExtendedButTreatedLikeSaurabh = DataFrame(envRobustAvg = colwise(mean, convert(DataFrame, transpose(convert(Array, envRobustnessDuringEvolutionExtendedButTreatedLikeSaurabh)))))
# # Plot results:
# # why is environmental robustness changing so much for the same population just different evolutionary trajectories?? -ML 4/3/19
# extendedPlot = plot(Vector(1:21),Array(envRobustnessDuringEvolutionExtended),title="Extended Model Env Robustness")
# extendedPlotTreatedLikeSaurabh = plot(Vector(1:21),Array(envRobustnessDuringEvolutionExtendedButTreatedLikeSaurabh),title="Extended Model Treated Saurabh Env Robustness")
# saurabhPlot = plot(Vector(1:21),Array(envRobustnessDuringEvolutionSaurabh),title="Saurabh Model Env Robustness")
# plot(extendedPlot,extendedPlotTreatedLikeSaurabh,saurabhPlot,layout=(3,1))
# extendedAvgPlot = plot(Vector(1:21),Array(averagesDuringEvolutionExtended),title="Extended Model Avg Env Robustness")
# extendedAvgPlotTreatedLikeSaurabh = plot(Vector(1:21),Array(averagesDuringEvolutionExtendedButTreatedLikeSaurabh),title="Extended Model Treated Saurabh Avg Env Robustness")
# saurabhAvgPlot = plot(Vector(1:21),Array(averagesDuringEvolutionSaurabh),title="Saurabh Model Avg Env Robustness")
# plot(extendedAvgPlot,extendedAvgPlotTreatedLikeSaurabh,saurabhAvgPlot,layout=(3,1))
# ################################################################################


##################################################################################
# Testing how many unstable perturbed individuals occur when testing for
# environmental robustness:
# Added these lines to environemntal robustness functions for our code and
# C++ code version:
# if !perturbed.stable
#     numStable = numStable + 1
# end
# print(" number unstable = ", numStable) at the end of function
saurabhPop = genpop()
founderGeneticRobustness = 0.
mutatedEnvMatrix = generateMutatedEnvsMatrix(saurabhPop.founder)
founderEnvRobustness = founderEnvRobust(saurabhPop.founder, mutatedEnvMatrix)
environmentalRobustnessSaurabhCcodeVersion(saurabhPop.individuals[1], mutatedEnvMatrix, founderEnvRobustness, saurabhPop.founder)
envRobustCcode = Float32(saurabhPop.individuals[1].envRobustness)
environmentalSensitivitySaurabhCcodeVersion(saurabhPop.individuals[1], mutatedEnvMatrix, saurabhPop.founder)
envSensCcode = Float32(saurabhPop.individuals[1].envRobustness)
environmentalRobustness(saurabhPop.individuals[1], mutatedEnvMatrix, founderEnvRobustness)
envRobust = Float32(saurabhPop.individuals[1].envRobustness)
develStateDiff = Float32(sum(abs.(saurabhPop.founder.develstate - saurabhPop.individuals[1].develstate))/G)
environmentalSensitivity(saurabhPop.individuals[1], mutatedEnvMatrix)
envSens = Float32(saurabhPop.individuals[1].envRobustness)
print(" Dif env robustness = ",envRobustCcode-envRobust)
print(" devel state diff btw founder & ind: ", develStateDiff)
print(" Dif env sens = ", envSensCcode-envSens)
print(" Dif btw robustness & devel difs = ", abs(envRobustCcode-envRobust)-develStateDiff)
for i = 1:200
    update(saurabhPop)
end
environmentalRobustnessSaurabhCcodeVersion(saurabhPop.individuals[1], mutatedEnvMatrix, founderEnvRobustness, saurabhPop.founder)
envRobustCcode = Float32(saurabhPop.individuals[1].envRobustness)
print(" env robustness C++ = ", envRobustCcode)
environmentalSensitivitySaurabhCcodeVersion(saurabhPop.individuals[1], mutatedEnvMatrix, saurabhPop.founder)
envSensCcode = Float32(saurabhPop.individuals[1].envRobustness)
environmentalRobustness(saurabhPop.individuals[1], mutatedEnvMatrix, founderEnvRobustness)
envRobust  = Float32(saurabhPop.individuals[1].envRobustness)
print(" env robustness = ", envRobust)
develStateDiff = Float32(sum(abs.(saurabhPop.founder.develstate - saurabhPop.individuals[1].develstate))/G)
environmentalSensitivity(saurabhPop.individuals[1], mutatedEnvMatrix)
envSens = Float32(saurabhPop.individuals[1].envRobustness)
print(" Dif env robustness = ",envRobustCcode-envRobust)
print(" devel state diff btw founder & ind: ", develStateDiff)
print(" Dif btw robustness & devel difs = ", abs(envRobustCcode-envRobust)-develStateDiff)
##################################################################################


@everywhere configfile = "constants.jl"
@everywhere indir = joinpath("..","input")
@everywhere constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
@everywhere include(constantsFile)
print(" SELSTR = ", SELSTR)
include("runMultipleIndependentTrials.jl")
SELSTR = 0.4
print(" SELSTR = ", SELSTR)
include("runMultipleIndependentTrials.jl")
SELSTR = 0.1
print(" SELSTR = ", SELSTR)
include("runMultipleIndependentTrials.jl")

POLYCOMBFLAG = true
SELSTR = 1.
print(" SELSTR = ", SELSTR)
include("runMultipleIndependentTrials.jl")
SELSTR = 0.4
print(" SELSTR = ", SELSTR)
include("runMultipleIndependentTrials.jl")
SELSTR = 0.1
print(" SELSTR = ", SELSTR)
include("runMultipleIndependentTrials.jl")
