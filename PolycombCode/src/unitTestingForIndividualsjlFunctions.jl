#############################################
# Unit Testing for individuals.jl functions
#############################################
using Test
using Distributions
using Dates
using DataFrames
#Pkg.build("GR") # if restarting julia
using Plots
using Combinatorics
using Printf
using Distributed


configfile = "constants.jl"
indir = joinpath("..","input")
constantsFile = joinpath(indir,configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

include("types.jl")
include("individuals.jl")
include("population.jl")
include("measure.jl")
include("utilities.jl")

################################################################################
# Test adding functionality for generating population when W matrix, and/or
# initstate, and/or optstate is given by user

### 1) When none given: (tested by Maryl on 3/25/19 and works)
const WMAT = Array{Float64,2}(undef,0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
popWithNone = genpop()
@test length(findall(x->x!=0,popWithNone.founder.network)) > 0
@test length(findall(x->x!=0,popWithNone.founder.optstate)) > 0
@test popWithNone.founder.develstate == popWithNone.founder.optstate
@test length(findall(x->x!=0,popWithNone.founder.initstate)) > 0
@test popWithNone.founder.stable == true

### 2) When initstate only given: (tested by Maryl on 3/25/19 and works)
const WMAT = Array{Float64,2}(undef,0,0)
const INP = convert(Vector{Float64},rand(0:1,G))
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
popWithInitstate = genpop()
@test INP == popWithInitstate.founder.initstate
@test WMAT != popWithInitstate.founder.network
@test length(findall(x->x!=0,popWithInitstate.founder.optstate)) > 0
@test popWithInitstate.founder.develstate == popWithInitstate.founder.optstate
@test popWithInitstate.founder.stable == true

### 3) When initstate & W matrix given: (tested by Maryl on 3/25/19 and works)
const WMAT = popWithInitstate.founder.network # will always work because this individual is already stable; if give random W matrix as WMAT then genpop() may not be able to get a stable founder --> error will return if so
const INP = popWithInitstate.founder.initstate
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
popWithInitAndWmat = genpop()
@test INP == popWithInitAndWmat.founder.initstate
@test WMAT == popWithInitAndWmat.founder.network
@test length(findall(x->x!=0,popWithInitAndWmat.founder.optstate)) > 0
@test popWithInitAndWmat.founder.develstate == popWithInitAndWmat.founder.optstate
@test popWithInitAndWmat.founder.stable == true

### 4) When W matrix only given: (tested by Maryl on 3/25/19 and works)
const WMAT = popWithNone.founder.network
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
popWithWmatOnly = genpop()
@test WMAT == popWithWmatOnly.founder.network
@test INP != popWithWmatOnly.founder.initstate
@test OPT != popWithWmatOnly.founder.optstate
@test length(findall(x->x!=0,popWithWmatOnly.founder.optstate)) > 0
@test popWithWmatOnly.founder.develstate == popWithWmatOnly.founder.optstate
@test popWithWmatOnly.founder.stable == true

#*** Fix these because doesn't need to find stable develstate that is same as
#*** the optstate or even close to the optstate

#**** These 3 cases should take forever to run because takes a long time to
# find an initstate and/or network (W mat) that leads to optstate and is stable
# --> but doesn't need to lead to optstate...so why taking forever? ** ML 4/3/19
### 5) When initstate & optstate given:
### 6) When optstate only given:
### 7) When W matrix & optstate given:

## Test with actual populations from Saurabh version to use to run extended version
SAURABVERSION = true
const WMAT = Array{Float64,2}(undef,0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
saurabhPop = genpop()

SAURABVERSION = false
const WMAT = saurabhPop.founder.network
const INP = saurabhPop.founder.initstate
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
extendedPop = genpop()
@test length(findall(x->x!=0,saurabhPop.founder.network - extendedPop.founder.network)) == 0 # should be zero because using newly implemented ability to generate population based on same W matrix (AKA network)
@test length(findall(x->x!=0,saurabhPop.founder.initstate - extendedPop.founder.initstate)) == 0 #should be zero
################################################################################


################################################################################
### Test iterateindfound() function: (tested by Maryl 3/26/19 and is good)
const WMAT = Array{Float64,2}(undef,0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
popTest = genpop()
popTestFounder = deepcopy(popTest.founder)
me = deepcopy(popTest.founder)
currstate = copy(me.initstate)
vectAD = zeros(Float64,MAXCONV)
expressionToPlot = DataFrame()
stateupdate = copy(me.initstate)
i = 1
for i=1:MAXCONV
    print(i)
    stateupdate = me.network*currstate
    stateupdateenv = me.EnvIndInter*me.EnvState
    stateupdate = stateupdate + stateupdateenv
    stateupdate = 1.0 ./ (1.0 .+ exp.(-SIGSTR.*stateupdate))
    tempdiff = currstate - stateupdate
    actualDist = sum(abs.(tempdiff))/G
    #sliding window - added 7/5/18 - SG
    vectAD[i] = actualDist
    if i >= TAU
        print(" >= TAU ")
        distBool = (vectAD[i] <= 10.0^(-6.0))
        stableBool = slidingWindow(i,vectAD)
        if distBool && stableBool
            print(" stable ")
            me.stable = true
            me.develstate = copy(stateupdate) #changed to copy() - SG 6/21/18
            me.pathlength = copy(i) #changed to copy() - SG 6/21/18
            expressionToPlot[Symbol(i)] = copy(stateupdate)
            break
        elseif i==MAXCONV
            print(" not stable ")
            me.stable = false
            me.develstate = copy(stateupdate) #changed to copy() - SG 6/21/18
            me.pathlength = copy(MAXCONV) #changed to copy() - SG 6/21/18
        end
    end
    currstate = copy(stateupdate)
    expressionToPlot[Symbol(i)] = copy(stateupdate)
end
@test popTestFounder.develstate == stateupdate # won't equal currstate because
# "break" call after reach stability happens before set currstate to stateupdate
@test popTestFounder.develstate == me.develstate
@test popTestFounder.pathlength == i
@test popTestFounder.pathlength == me.pathlength
plot(Vector(1:size(expressionToPlot,2)),transpose(convert(Matrix,expressionToPlot)),legend=false)
################################################################################
### Test slidingWindow() function: (tested by Maryl 3/26/19 and is good)
stds = []
for t = 1:me.pathlength # me.pathlength is t where found stable by sliding window
    x = vectAD
    if t >= TAU
        dists = sum(x[(t-(TAU-1)):t]) #add up each of the distances between measurements
        push!(stds, dists/TAU) #find std by divding by TAU
    end
end
plot(stds)
plot!([EPSI], seriestype="hline")
################################################################################

################################################################################
# Added ability and testing taking measurements based on user input:
vecMeasurementsToTake = ["connectivity","envRobustness"]#ARGS[3]
data = Dict{String, DataFrame}()
dataTestMeasureTypes = Dict{String, DataFrame}()
for i=1:length(vecMeasurementsToTake)
    data[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
    dataTestMeasureTypes[string(vecMeasurementsToTake[i],"DuringEvolution")] = DataFrame()
end
const WMAT = Array{Float64,2}(undef,0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
include("measure.jl")
#for indTrial=1:2
    popTest = genpop()
    meas = genmeasure()
    measureTypes = Dict("time"=>meas.time, "fitness"=>meas.fitness, "fitnessstd"=>meas.fitnessstd,
                    "robustness"=>meas.robustness, "robustnessstd"=>meas.robustnessstd, "envRobustness"=>meas.envRobustness,
                    "envRobustnessStd"=>meas.envRobustnessStd, "pathlength"=>meas.pathlength, "pathlengthstd"=>meas.pathlengthstd,
                    "indtypes"=>meas.indtypes, "inittypes"=>meas.inittypes, "develtypes"=>meas.develtypes, "pcgVecTypes"=>meas.pcgVecTypes,
                    "pcgStateTypes"=>meas.pcgStateTypes, "opttypes"=>meas.opttypes, "connectivity"=>meas.connectivity)
    mutatedEnvMatrix = generateMutatedEnvsMatrix(popTest.founder)
    founderGeneticRobust(popTest.founder)
    founderGeneticRobustness = popTest.founder.robustness
    founderEnvRobust(popTest.founder)
    founderEnvRobustness = popTest.founder.envRobustness
    measure(popTest, meas, 1, 1, mutatedEnvMatrix, founderGeneticRobustness,
        founderEnvRobustness, vecMeasurementsToTake, measureTypes)
    update(popTest)
    measure(popTest, meas, 2, 2, mutatedEnvMatrix, founderGeneticRobustness,
        founderEnvRobustness, vecMeasurementsToTake, measureTypes)
    # Unit testing to make sure meas.measurement equals measureTypes[measurement],
    # so when set for measureTypes[measurement] that meas.measurement is also updated
    for i=1:length(vecMeasurementsToTake)
        print("\n", string(vecMeasurementsToTake[i])," ", @test  getfield(meas,Symbol(vecMeasurementsToTake[i])) == measureTypes[vecMeasurementsToTake[i]])
        print("\n", string(vecMeasurementsToTake[i])," equals ", getfield(meas,Symbol(vecMeasurementsToTake[i]))[1:2]) # [1:2] bc only doing update twice for this test
    end
    for i=1:length(vecMeasurementsToTake)
        data[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = getfield(meas,Symbol(vecMeasurementsToTake[i]))
        dataTestMeasureTypes[string(vecMeasurementsToTake[i],"DuringEvolution")][Symbol(string("trial",indTrial))] = measureTypes[vecMeasurementsToTake[i]]
    end
    # Making sure measurement values don't change:
    for i=1:length(vecMeasurementsToTake)
        print("\n", string(vecMeasurementsToTake[i])," equals after save to dataframe ", getfield(meas,Symbol(vecMeasurementsToTake[i]))[1:2]) # [1:2] bc only doing update twice for this test
    end
#end
for i=1:length(vecMeasurementsToTake)
    data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average] = colwise(mean, convert(DataFrame, transpose(convert(Matrix, data[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
    data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Sem] = colwise(standardErrorOfMean, convert(DataFrame, transpose(convert(Matrix, data[string(vecMeasurementsToTake[i],"DuringEvolution")]))))
end
averagesDuringEvolution = DataFrame()
for i=1:length(vecMeasurementsToTake)
    averagesDuringEvolution[Symbol(string(vecMeasurementsToTake[i],"Avg"))] = data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Average]
    averagesDuringEvolution[Symbol(string(vecMeasurementsToTake[i],"Sem"))] = data[string(vecMeasurementsToTake[i],"DuringEvolution")][:Sem]
end


################################################################################
# Test generateMutatedEnvsMatrix() function
################################################################################
SAURABVERSION = true
const WMAT = Array{Float64,2}(undef,0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
saurabhPop = genpop()
mutatedInitStates = collect(combinations(collect(1:length(saurabhPop.founder.initstate)), convert(Int64, ALTERENV*G)))
mutatedInitStateIndx = sample(1:length(mutatedInitStates), ENVCHANGES, replace = false)
## Test way find random entries to mutate:
differences = Array(1:ENVCHANGES)
for i = 1:ENVCHANGES
    differences[i] = diff(mutatedInitStates[mutatedInitStateIndx[i]])[1]
end
differencesAlter = Array(1:ENVCHANGES)
mutatedInitStatesAlter = zeros(Int64,(ENVCHANGES,2))
for i = 1:ENVCHANGES
    mutatedInitStatesAlter[i,1] = rand(1:length(saurabhPop.founder.initstate))
    secondState = rand(1:length(saurabhPop.founder.initstate))
    while true
        mutatedInitStatesAlter[i,2] = secondState
        secondState != mutatedInitStatesAlter[i,1] && break
        secondState = rand(1:length(saurabhPop.founder.initstate))
    end
    differencesAlter[i] = abs(diff(mutatedInitStatesAlter[i,:])[1])
end
plot(histogram(differences, title="differences"),histogram(differencesAlter, title="differencesAlter"),layout=(2,1))
## Way calculate random entries to choose to mutate is fine as originally was -ML 4/4/19


###############################################################################
# Test calculating founder's environmental robustness:
###############################################################################
local tempdiff
dist = 0.0
if SAURABVERSION
    for i=1:ENVCHANGES
        perturbed = deepcopy(me)
        perturbed.initstate = copy(mutatedEnvMatrix[:,i])
        iterateind(perturbed) # use this function to test if stability is reached
        tempdiff = me.develstate - perturbed.develstate
        dist += sum(abs.(tempdiff))/(G)
    end
else
    for i=1:ENVCHANGES
        perturbed = deepcopy(me)
        perturbed.EnvState = copy(mutatedEnvMatrix[:,i])
        iterateind(perturbed) # use this function to test if stability is reached
        tempdiff = me.develstate - perturbed.develstate
        dist += sum(abs.(tempdiff))/(G)
    end
end
me.envRobustness = dist/ENVCHANGES


################################################################################
# Test environmental robustness behavior before evolution
################################################################################
### With Saurabh version (tests pass -ML 4/4/19)
SAURABVERSION = true
const WMAT = Array{Float64,2}(undef,0,0)
const INP = (Float64)[]
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
saurabhPop = genpop()

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
if "robustness" in vecMeasurementsToTake
    founderGeneticRobustness = founderGeneticRobust(saurabhPop.founder)
else
    founderGeneticRobustness = 0.
end
if "envRobustness" in vecMeasurementsToTake
    mutatedEnvMatrix = generateMutatedEnvsMatrix(saurabhPop.founder)
    founderEnvRobustness = founderEnvRobust(saurabhPop.founder)
else
    founderEnvRobustness = 0.
    mutatedEnvMatrix = Array{Float64,2}[]
end
t = 1 # at first generation
measure(saurabhPop, meas, t, measnum, mutatedEnvMatrix, founderGeneticRobustness, founderEnvRobustness, vecMeasurementsToTake, measureTypes)
@test meas.envRobustness[1] == 0 #average of population robustness is 0 because haven't evolved yet
@test environmentalRobustness(saurabhPop.founder,mutatedEnvMatrix,founderEnvRobustness) == 0 #robustness of founder should be zero b/c subtracting from founder

### With extended version: (tests pass -ML 4/4/19)
SAURABVERSION = false
const WMAT = saurabhPop.founder.network
const INP = saurabhPop.founder.initstate
const OPT = (Float64)[]
include("types.jl")
include("individuals.jl")
include("population.jl")
extendedPop = genpop()
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
if "robustness" in vecMeasurementsToTake
    founderGeneticRobustness = founderGeneticRobust(extendedPop.founder)
else
    founderGeneticRobustness = 0.
end
if "envRobustness" in vecMeasurementsToTake
    mutatedEnvMatrix = generateMutatedEnvsMatrix(extendedPop.founder)
    founderEnvRobustness = founderEnvRobust(extendedPop.founder)
else
    founderEnvRobustness = 0.
    mutatedEnvMatrix = Array{Float64,2}[]
end
t = 1 # at first generation
measure(extendedPop, meas, t, measnum, mutatedEnvMatrix, founderGeneticRobustness, founderEnvRobustness, vecMeasurementsToTake, measureTypes)
@test meas.envRobustness[1] == 0 #average of population robustness is 0 because haven't evolved yet
environmentalRobustness(extendedPop.founder,mutatedEnvMatrix,founderEnvRobustness)
@test extendedPop.founder.envRobustness == 0 #robustness of founder should be zero b/c subtracting from founder
################################################################################



################################################################################
# Test mutate() function:
################################################################################
connectivityOfW = C*G^2
popTest = genpop()
nzindx = findall(x->x!=0,popTest.founder.network)
@test popTest.founder.network[nzindx[1]] != 0 # make sure index of findall is working properly
cnum = length(nzindx)
print("\n", "Desired connectivity: ", connectivityOfW," vs actual connectivity: ", cnum)

mutflag = zeros(10000)
for i = 1:10000
    mutflag[i] = mutate(popTest.founder)[1]
end
mutationRate = length(findall(x->x!=0,mutflag))/length(mutflag)
print("\n", "Desired mutation rate: ", MUTRATE," vs actual rate: ", mutationRate) #check to make sure close to 10% (0.1)

pop = genpop()
popTest = deepcopy(pop)
mutflag = 0
mutatedElement = 0
while mutflag == 0
    mutflag, mutatedElement = mutate(popTest.founder)
end
@test popTest.founder.network != pop.founder.network #make sure network has been changed with mutation
@test popTest.founder.network[mutatedElement] != pop.founder.network[mutatedElement] # make sure mutatedElement actually has been element to be changed

# Make sure mutations are follow normal distribution:
mutationChange = zeros(1000)
mutationValue = zeros(1000)
for i = 1:1000
    pop = genpop()
    popTest = deepcopy(pop)
    mutflag = 0
    while mutflag == 0
        mutflag, mutatedElement = mutate(popTest.founder)
    end
    mutationChange[i] = popTest.founder.network[mutatedElement] - pop.founder.network[mutatedElement]
    mutationValue[i] = popTest.founder.network[mutatedElement]
end
p = plot(sort(mutationChange))
plot!(p,sort(randn(1000))) # mutation changes do NOT follow a normal distribution, but they don't necessarily have to, right? -ML 5/2/19
mean(mutationChange) # mean is not that close to 0...
std(mutationChange) # std is not that close to 1...

plotVal = plot(sort(mutationValue))
plot!(plotVal,sort(randn(1000)))
mean(mutationValue)
std(mutationValue)
################################################################################


################################################################################
# Test connectivity of EnvIndInter matrix (W x E)
################################################################################
connectivityOfEnvW = (CENMAT*CENGENE)*(G*ENVS) # currently (5/2/19) code is written such that connectivity of EnvIndInter is C*C and not C like W matrix.... -ML --> Aviv says this is OK and connectivityOfEnvW could even be smaller than 0.16...
popTest = genpop()
nzindx = findall(x->x!=0,popTest.founder.EnvIndInter)
# make sure no environmental interactions equal 0 that were changed:
@test length(findall(y->y==0,map(x->popTest.founder.EnvIndInter[nzindx[x]],1:length(nzindx)))) == 0
# look at connectivity of environmental interaction matrix of actual vs desired:
cnum = length(nzindx)
print("\n", "Desired connectivity: ", connectivityOfEnvW," vs actual connectivity: ", cnum)
## look at what's happening with genes and environments in total separately, so columns and rows
# Environments affecting genes (rows):
proportionOfEnvsAffectingGene=map(y->length(findall(x->x!=0,popTest.founder.EnvIndInter[y,:]))/G,1:ENVS)
# Proportion of how many environments affect at least 1 gene: (should be around CENMAT which is 0.4)
length(findall(x->x!=0,proportionOfEnvsAffectingGene))/G
# Average number of genes a given environment affects: (should be around CENGENE which is 0.4 because environment must pass CENMAT in order to affect any genes)
mean(proportionOfEnvsAffectingGene[findall(x->x!=0,proportionOfEnvsAffectingGene)])
# Genes affecting environments (columns):
proportionOfGenesAffectingEnv=map(y->length(findall(x->x!=0,popTest.founder.EnvIndInter[:,y]))/ENVS,1:G)
# Proportion of how many genes are affected by at least 1 enivonrment:
length(findall(x->x!=0,proportionOfGenesAffectingEnv))/ENVS # = 1.00 because in ENVS times randn() will be less than CENGENE at least once
# Average number of environments that affect a given gene: (should be around 0.16 (= CENMAT * CENGENE))
mean(proportionOfGenesAffectingEnv[findall(x->x!=0,proportionOfGenesAffectingEnv)])

### TEST by changing CENMAT and CENGENE: ###
CENMAT = 0.2
CENGENE = 0.3
include("individuals.jl")
connectivityOfEnvW = (CENMAT*CENGENE)*(G*ENVS) # currently (5/2/19) code is written such that connectivity of EnvIndInter is C*C and not C like W matrix.... -ML --> Aviv says this is OK and connectivityOfEnvW could even be smaller than 0.16...
popTest = genpop()
nzindx = findall(x->x!=0,popTest.founder.EnvIndInter)
# make sure no environmental interactions equal 0 that were changed:
@test length(findall(y->y==0,map(x->popTest.founder.EnvIndInter[nzindx[x]],1:length(nzindx)))) == 0
# look at connectivity of environmental interaction matrix of actual vs desired:
cnum = length(nzindx)
print("\n", "Desired connectivity: ", connectivityOfEnvW," vs actual connectivity: ", cnum)
## look at what's happening with genes and environments in total separately, so columns and rows
# Environments affecting genes (rows):
proportionOfEnvsAffectingGene=map(y->length(findall(x->x!=0,popTest.founder.EnvIndInter[y,:]))/G,1:ENVS)
# Proportion of how many environments affect at least 1 gene: (should be around CENMAT which is 0.2)
length(findall(x->x!=0,proportionOfEnvsAffectingGene))/G
# Average number of genes a given environment affects: (should be around CENGENE which is 0.3 because environment must pass CENMAT in order to affect any genes)
mean(proportionOfEnvsAffectingGene[findall(x->x!=0,proportionOfEnvsAffectingGene)])
# Genes affecting environments (columns):
proportionOfGenesAffectingEnv=map(y->length(findall(x->x!=0,popTest.founder.EnvIndInter[:,y]))/ENVS,1:G)
# Proportion of how many genes are affected by at least 1 enivonrment:
length(findall(x->x!=0,proportionOfGenesAffectingEnv))/ENVS # = 1.00 because in ENVS times randn() will be less than CENGENE at least once
# Average number of environments that affect a given gene: (should be around 0.06 (= CENMAT * CENGENE = 0.2 * 0.3))
mean(proportionOfGenesAffectingEnv[findall(x->x!=0,proportionOfGenesAffectingEnv)])
# --> tests and assumptions are all correct! 5/23/19
################################################################################


################################################################################
# Test mutateEnvIndInter() function
################################################################################
envmutflag = zeros(1000)
for i = 1:1000
    envmutflag[i] = mutateEnvIndInter(popTest.founder)[1]
end
envMutationRate = length(findall(x->x!=0,envmutflag))/length(envmutflag)
print("\n", "Desired env mutation rate: ", ENVMUTRATE," vs actual: ", envMutationRate) #check to make sure close to 10% (0.1)

pop = genpop()
popTest = deepcopy(pop)
mutflag = 0
mutatedElement = 0
while mutflag == 0
    mutflag, mutatedElement = mutateEnvIndInter(popTest.founder)
end
@test popTest.founder.EnvIndInter != pop.founder.EnvIndInter #make sure network has been changed with mutation
@test popTest.founder.EnvIndInter[mutatedElement] != pop.founder.EnvIndInter[mutatedElement] # make sure mutatedElement actually has been element to be changed

# fix so tests mutateEnvIndInter instead of mutate (did this above for mutate function)
# Make sure mutations are follow normal distribution:
mutationChange = zeros(1000)
mutationValue = zeros(1000)
for i = 1:1000
    pop = genpop()
    popTest = deepcopy(pop)
    mutflag = 0
    mutatedElement = 0
    while mutflag == 0
        mutflag, mutatedElement = mutateEnvIndInter(popTest.founder)
    end
    mutationChange[i] = popTest.founder.EnvIndInter[mutatedElement] - pop.founder.EnvIndInter[mutatedElement]
    mutationValue[i] = popTest.founder.EnvIndInter[mutatedElement]
end
p = plot(sort(mutationChange))
plot!(p,sort(randn(1000))) # mutation changes do follow a normal distribution
mean(mutationChange) # mean is close to 0...
std(mutationChange) # std is close to 1...

plotVal = plot(sort(mutationValue))
plot!(plotVal,sort(randn(100)))
mean(mutationValue)
std(mutationValue)
################################################################################


################################################################################
# Test polymut() function
################################################################################
# test that polycombvec doesn't change during polymut() function except where
# supposed to be mutated in function
include("individuals.jl")
pop = genpop()
popTest = deepcopy(pop)
polymutIndices = []
runsToPolyMutation = 0
while length(polymutIndices) == 0 # run polymut() function until mutation occurs
    # and keep track of how many times run before get polycomb mutation
    polymutVec = polymut(popTest.founder)
    polymutIndices = findall(x->x!=0,polymutVec)
    runsToPolyMutation = runsToPolyMutation + 1
end
print(runsToPolyMutation)
@test map(x->popTest.founder.polycombvec[polymutIndices[x]],1:length(polymutIndices)) != map(x->pop.founder.polycombvec[polymutIndices[x]],1:length(polymutIndices))
################################################################################


################################################################################
# Test fitnesseval() function
################################################################################
pop = genpop()
popTest = deepcopy(pop)
statediff = popTest.individuals[1].optstate - popTest.individuals[1].develstate
distance= sum(abs.(statediff))/(G) # ranges from 0 to infinity
fitnessOfIndividual = exp.(-(distance/SELSTR))
# Before evolution individual's optstate and individual's develstate should be the same:
@test distance == 0
@test fitnessOfIndividual == 1

for i = 1:50
    update(popTest)
end
statediffAfterEvolution = popTest.individuals[1].optstate - popTest.individuals[1].develstate
distanceAfterEvolution = sum(abs.(statediffAfterEvolution))/(G) # ranges from 0 to infinity
fitnessOfIndividualAFterEvolution = exp.(-(distanceAfterEvolution/SELSTR))
# now individual's optstate and individual's develstate should be different
# so fitness is no longer 1:
@test distanceAfterEvolution != 0
@test fitnessOfIndividualAFterEvolution != 1
################################################################################


################################################################################
# Test founderGeneticRobust() function
################################################################################
include("individuals.jl")
pop = genpop()
popTest = deepcopy(pop)

repeatedTimes = 1000
stablePerturbations = zeros(repeatedTimes)
founderGeneticRobustness = zeros(repeatedTimes)
ROBIT = 2500
include("individuals.jl")
for i = 1:repeatedTimes
    stablePerturbations[i] = founderGeneticRobust(popTest.founder)
    founderGeneticRobustness[i] = popTest.founder.robustness
end
mean(stablePerturbations)
std(stablePerturbations)
mean(founderGeneticRobustness)
std(founderGeneticRobustness)
# Does genetic robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot(stablePerturbations, founderGeneticRobustness, seriestype=:scatter) # seems like more stable
# perturbations the lower the genetic robustness is -ML 5/29/19
# **Should we only measure robustness for the those perturbations that produce
# **stable individuals?? -ML 5/29/19 --> Aviv says no, that if anything we
# should do more perturbations (so increase ROBIT; tested with 100, 1000, 2500
# in code above and saved results in word doc called:
# founderGeneticRobustnessUnitTestGraphs.docx) -6/3/19
################################################################################


################################################################################
# Test geneticRobustness() function
################################################################################
pop = genpop()
popTest = deepcopy(pop)

ROBIT = 100
include("individuals.jl")
# test 10,000 perturbations vs 100 perturbations when calculating founder's genetic robustness:
include("individuals.jl")
founderStablePerturbations = zeros(1000)
founderGeneticRobustness = zeros(1000)
for j = 1:1000
    founderStablePerturbations[j] = founderGeneticRobust(pop.founder)
    founderGeneticRobustness[j] = pop.founder.robustness
end
plot(founderGeneticRobustness, seriestype=:scatter, title="10,000 perturbations")
plot(founderStablePerturbations, founderGeneticRobustness, seriestype=:scatter,title="10,000 perturbations")
# test robustness of individual's in population and # of stable perturbations:
founderGeneticRobust(pop.founder)
print("founder robustness = ", pop.founder.robustness)
stablePerturbationsBeforeEvolution = map(x->geneticRobustness(pop.individuals[x], pop.founder.robustness), 1:N)
geneticRobustnessBeforeEvolution = map(x->pop.individuals[x].robustness, 1:N)
#plot(stablePerturbationsBeforeEvolution, geneticRobustnessBeforeEvolution,seriestype=:scatter, title="all individuals before evolution")
meanRobustnessBeforeEvolution = mean(geneticRobustnessBeforeEvolution)
print(" mean = ", meanRobustnessBeforeEvolution)
print(" median = ", median(geneticRobustnessBeforeEvolution)) # perhaps should do median so not as affected by outliers?
print(" mean stable perturbations = ", mean(stablePerturbationsBeforeEvolution)) # If rerun founderGeneticRobust(popTest.founder)
# then can change the mean & median, but if do not rerun that but rerun all after that and above this, then get similar results....
# is this ok? It happens before and after evolution (well 100 generations). -ML 6/10/19
for j = 1:100 # undergo 100 generations of evolution
    update(pop)
end
robustnessDifferentFromFounder = length(findall(x->x!=0,map(x->geneticRobustness(pop.individuals[x], pop.founder.robustness), 1:N)))
print(robustnessDifferentFromFounder)
@test robustnessDifferentFromFounder > 0 #robustness of individuals should not be zero anymore because should be different from the founder
individualNum = 500
meanStablePerturbations = zeros(individualNum)
meanGeneticRobustness = zeros(individualNum)
repeatedTimes = 100
stablePerturbations = zeros(repeatedTimes)
individualGeneticRobustness = zeros(repeatedTimes)
for k = 1:individualNum
    for i = 1:repeatedTimes
        stablePerturbations[i] = geneticRobustness(pop.individuals[individualNum], pop.founder.robustness)
        individualGeneticRobustness[i] = pop.individuals[individualNum].robustness
    end
    meanStablePerturbations[k] = mean(stablePerturbations)
    std(stablePerturbations)
    meanGeneticRobustness[k] = mean(individualGeneticRobustness)
    std(individualGeneticRobustness)
end
meanOfIndividuals = mean(meanGeneticRobustness)
print("mean after 100 generations = ", meanOfIndividuals)
print(" median after 100 generations = ", median(meanGeneticRobustness))
print(" mean stable perturbations = ", mean(meanStablePerturbations))
# Does genetic robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot(stablePerturbations, individualGeneticRobustness, seriestype=:scatter)
# For all individuals: still correlates; basically same test as with founder but after running evolution for 100 generations:
plot(meanStablePerturbations, meanGeneticRobustness, seriestype=:scatter, title="mean of individuals")
################################################################################


################################################################################
# Test environmental robustness:
################################################################################
pop = genpop()
popTest = deepcopy(pop)

ENVCHANGES = 1225 # 1225 is all the possible combinations of changing 2 of the 50 environmental states (50 choose 2)
include("individuals.jl")
# test 10,000 perturbations vs 100 perturbations when calculating founder's environmental robustness:
include("individuals.jl")
founderStablePerturbations = zeros(100)
founderEnvRobustness = zeros(100)
for j = 1:100 # test with different randomly choosen perturbed environments and
    # see how much varies depending on number of perturbed environments tested:
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    founderStablePerturbations[j] = founderEnvRobust(pop.founder) # *** add measuring of stable perturbations to env robustness functions!! -ML 6/11/19
    founderEnvRobustness[j] = pop.founder.envRobustness
end
print("mean env robustness = ", mean(founderEnvRobustness), " std = ", std(founderEnvRobustness))
plot(founderEnvRobustness, seriestype=:scatter, title="1,000 perturbations; Env Robustness", xaxis="# founders tested", yaxis="Env Robustness")
plot(founderStablePerturbations, founderEnvRobustness, seriestype=:scatter,title="1,000 perturbations; Env Robustness",xaxis="# stable perturbations",yaxis="Env Robustness")
# test robustness of individual's in population and # of stable perturbations:
founderEnvRobust(pop.founder)
print("founder env robustness = ", pop.founder.envRobustness)
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbationsBeforeEvolution = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
envRobustnessBeforeEvolution = map(x->pop.individuals[x].envRobustness, 1:N)
#plot(stablePerturbationsBeforeEvolution, envRobustnessBeforeEvolution,seriestype=:scatter, title="all individuals before evolution")
meanRobustnessBeforeEvolution = mean(envRobustnessBeforeEvolution)
print(" mean before evolution for population = ", meanRobustnessBeforeEvolution)
print(" median before evolution for population = ", median(envRobustnessBeforeEvolution)) # perhaps should do median so not as affected by outliers?
print(" mean stable perturbations before evolution for population = ", mean(stablePerturbationsBeforeEvolution)) # If rerun founderEnvRobust(popTest.founder)
# then can change the mean & median, but if do not rerun that but rerun all after that and above this, then get similar results....
# is this ok? It happens before and after evolution (well 100 generations). -ML 6/10/19
for j = 1:100 # undergo 100 generations of evolution
    update(pop)
end
robustnessDifferentFromFounder = length(findall(x->x!=0,map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)))
print(robustnessDifferentFromFounder)
@test robustnessDifferentFromFounder > 0 #robustness of individuals should not be zero anymore because should be different from the founder

ENVCHANGES = 100
include("individuals.jl")
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbations = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
individualEnvRobustness = map(x->pop.individuals[x].envRobustness, 1:N)
print("mean env robustness after 100 generations with ENVCHANGES 100 = ", mean(individualEnvRobustness), " and std = ", std(individualEnvRobustness))
print(" median after 100 generations = ", median(individualEnvRobustness))
print(" mean stable perturbations with ENVCHANGES 100 = ", mean(stablePerturbations), " and std = ", std(stablePerturbations))
# Does environmental robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot100envchanges = plot(stablePerturbations, individualEnvRobustness, seriestype=:scatter, title="env. robustness of individuals ENVCHANGES=100 & gens=100")

ENVCHANGES = 700
include("individuals.jl")
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbations = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
individualEnvRobustness = map(x->pop.individuals[x].envRobustness, 1:N)
print("mean env robustness after 100 generations with ENVCHANGES 700 = ", mean(individualEnvRobustness), " and std = ", std(individualEnvRobustness))
print(" median after 100 generations = ", median(individualEnvRobustness))
print(" mean stable perturbations with ENVCHANGES 700 = ", mean(stablePerturbations), " and std = ", std(stablePerturbations))
# Does environmental robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot700envchanges = plot(stablePerturbations, individualEnvRobustness, seriestype=:scatter, title="env. robustness of individuals ENVCHANGES=700 & gens=100")

ENVCHANGES = 1224
include("individuals.jl")
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbations = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
individualEnvRobustness = map(x->pop.individuals[x].envRobustness, 1:N)
print("mean env robustness after 100 generations with ENVCHANGES 1224 = ", mean(individualEnvRobustness), " and std = ", std(individualEnvRobustness))
print(" median after 100 generations = ", median(individualEnvRobustness))
print(" mean stable perturbations with ENVCHANGES 1224 = ", mean(stablePerturbations), " and std = ", std(stablePerturbations))
# Does environmental robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot1224envchanges = plot(stablePerturbations, individualEnvRobustness, seriestype=:scatter, title="env. robustness of individuals ENVCHANGES=1224 & gens=100")

# Now do 500 generations:
for j = 1:500 # undergo 100 generations of evolution
    update(pop)
end
ENVCHANGES = 100
include("individuals.jl")
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbations = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
individualEnvRobustness = map(x->pop.individuals[x].envRobustness, 1:N)
print("mean env robustness after 600 generations with ENVCHANGES 100 = ", mean(individualEnvRobustness), " and std = ", std(individualEnvRobustness))
print(" median after 600 generations = ", median(individualEnvRobustness))
print(" mean stable perturbations with ENVCHANGES 100 = ", mean(stablePerturbations), " and std = ", std(stablePerturbations))
# Does environmental robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot100envchanges600gens = plot(stablePerturbations, individualEnvRobustness, seriestype=:scatter, title="env. robustness of individuals ENVCHANGES=100 & gens=600")

ENVCHANGES = 700
include("individuals.jl")
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbations = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
individualEnvRobustness = map(x->pop.individuals[x].envRobustness, 1:N)
print("mean env robustness after 600 generations with ENVCHANGES 700 = ", mean(individualEnvRobustness), " and std = ", std(individualEnvRobustness))
print(" median after 600 generations = ", median(individualEnvRobustness))
print(" mean stable perturbations with ENVCHANGES 700 = ", mean(stablePerturbations), " and std = ", std(stablePerturbations))
# Does environmental robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot700envchanges600gens = plot(stablePerturbations, individualEnvRobustness, seriestype=:scatter, title="env. robustness of individuals ENVCHANGES=700 & gens=600")

ENVCHANGES = 1224
include("individuals.jl")
mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
stablePerturbations = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
individualEnvRobustness = map(x->pop.individuals[x].envRobustness, 1:N)
print("mean env robustness after 600 generations with ENVCHANGES 1224 = ", mean(individualEnvRobustness), " and std = ", std(individualEnvRobustness))
print(" median after 600 generations = ", median(individualEnvRobustness))
print(" mean stable perturbations with ENVCHANGES 1224 = ", mean(stablePerturbations), " and std = ", std(stablePerturbations))
# Does environmental robustness value correlate with # number of stable perturbations? --> YES -ML 6/3/19
plot1224envchanges600gens = plot(stablePerturbations, individualEnvRobustness, seriestype=:scatter, title="env. robustness of individuals ENVCHANGES=1224 & gens=600")

##########################################################
### Test 100 vs 1224 ENVCHANGES throughout evolution:
# Initialize number of trials to run for different populations to compare:
numTests = 2
# Initialize dataframes to store results:
resultsDataFrame100changes = DataFrame(meanChangesBeforeEvol = Float64[], meanChanges100gens = Float64[], meanChanges600gens = Float64[])
resultsDataFrame1224changes = DataFrame(meanChangesBeforeEvol = Float64[], meanChanges100gens = Float64[], meanChanges600gens = Float64[])
stdResultsDataFrame100changes = DataFrame(stdChangesBeforeEvol = Float64[], stdChanges100gens = Float64[], stdChanges600gens = Float64[])
stdResultsDataFrame1224changes = DataFrame(stdChangesBeforeEvol = Float64[], stdChanges100gens = Float64[], stdChanges600gens = Float64[])
stablePerturbationsResults100changes = DataFrame(stablePerturbationsBeforeEvol = Float64[], stablePerturbations100gens = Float64[], stablePerturbations600gens = Float64[])
stablePerturbationsResults1224changes = DataFrame(stablePerturbationsBeforeEvol = Float64[], stablePerturbations100gens = Float64[], stablePerturbations600gens = Float64[])
#
for i = 1:numTests
    pop = genpop()
    # before evolution:
    founderEnvRobust(pop.founder)
    ENVCHANGES = 100
    include("individuals.jl")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbations100changesBeforeEvol = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustness100changesBeforeEvol = map(x->pop.individuals[x].envRobustness, 1:N)

    ENVCHANGES = 1224
    include("individuals.jl")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbations1224changesBeforeEvol = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustness1224changesBeforeEvol = map(x->pop.individuals[x].envRobustness, 1:N)

    # after 100 generations of evolution:
    for j = 1:100
        update(pop)
    end
    ENVCHANGES = 100
    include("individuals.jl")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbations100changes100gens = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustness100changes100gens = map(x->pop.individuals[x].envRobustness, 1:N)

    ENVCHANGES = 1224
    include("individuals.jl")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbations1224changes100gens = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustness1224changes100gens = map(x->pop.individuals[x].envRobustness, 1:N)

    # after 600 generations of evolution:
    for j = 1:500
        update(pop)
    end
    ENVCHANGES = 100
    include("individuals.jl")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbations100changes600gens = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustness100changes600gens = map(x->pop.individuals[x].envRobustness, 1:N)

    ENVCHANGES = 1224
    include("individuals.jl")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbations1224changes600gens = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustness1224changes600gens = map(x->pop.individuals[x].envRobustness, 1:N)
    # Add all the results to dataframes:
    push!(resultsDataFrame100changes, [mean(individualEnvRobustness100changesBeforeEvol),mean(individualEnvRobustness100changes100gens),mean(individualEnvRobustness100changes600gens)])
    push!(resultsDataFrame1224changes, [mean(individualEnvRobustness1224changesBeforeEvol),mean(individualEnvRobustness1224changes100gens),mean(individualEnvRobustness1224changes600gens)])
    push!(stdResultsDataFrame100changes, [std(individualEnvRobustness100changesBeforeEvol),std(individualEnvRobustness100changes100gens),std(individualEnvRobustness100changes600gens)])
    push!(stdResultsDataFrame1224changes, [std(individualEnvRobustness1224changesBeforeEvol),std(individualEnvRobustness1224changes100gens),std(individualEnvRobustness1224changes600gens)])
    push!(stablePerturbationsResults100changes, [mean(stablePerturbations100changesBeforeEvol),mean(stablePerturbations100changes100gens),mean(stablePerturbations100changes600gens)])
    push!(stablePerturbationsResults1224changes, [mean(stablePerturbations1224changesBeforeEvol),mean(stablePerturbations1224changes100gens),mean(stablePerturbations1224changes600gens)])
end
# Compare throughout evolution when have 100 environmental changes and 1224 environmental changes separately:
using StatsPlots
plot100changes = @df resultsDataFrame100changes plot([:meanChangesBeforeEvol,:meanChanges100gens,:meanChanges600gens],title="100 envChanges comparisons",label=["Before Evol","100 gens","600 gens"],lw=2,ylabel="env robustness")#,yerror=stdResultsDataFrame100changes[:stdChangesBeforeEvol])
plot1224changes = @df resultsDataFrame1224changes plot([:meanChangesBeforeEvol,:meanChanges100gens,:meanChanges600gens],title="1224 envChanges comparisons",label=["Before Evol","100 gens","600 gens"],lw=2,xlabel="different population trials",ylabel="env robustness")#,yerror=stdResultsDataFrame1224changes)
plot(plot100changes, plot1224changes, layout=(2,1))

# Compare 100 vs 1224 environemtnal changes at each measured step in evolution:
plot100vs1224changesBeforeEvol = plot([resultsDataFrame100changes[:meanChangesBeforeEvol],resultsDataFrame1224changes[:meanChangesBeforeEvol]],title="EnvChanges comparison before evolution",label=["100 changes","1224 changes"],lw=2,ylabel="env robustness")
plot100vs1224changes100gens = plot([resultsDataFrame100changes[:meanChanges100gens],resultsDataFrame1224changes[:meanChanges100gens]],title="EnvChanges comparison after 100 gens",label=["100 changes","1224 changes"],lw=2,ylabel="env robustness")
plot100vs1224changes600gens = plot([resultsDataFrame100changes[:meanChanges600gens],resultsDataFrame1224changes[:meanChanges600gens]],title="EnvChanges comparison after 600 gens",label=["100 changes","1224 changes"],lw=2,xlabel="different population trials",ylabel="env robustness")
plot(plot100vs1224changesBeforeEvol,plot100vs1224changes100gens,plot100vs1224changes600gens, layout=(3,1))

# Print out means of the trials for each case:
print("\nmean env Robustness before evolution when 100 env changes = ", round(mean(resultsDataFrame100changes[:meanChangesBeforeEvol]); digits=7))
print("\nmean env Robustness after 100 gens when 100 env changes = ", round(mean(resultsDataFrame100changes[:meanChanges100gens]); digits=7))
print("\nmean env Robustness after 600 gens when 100 env changes = ", round(mean(resultsDataFrame100changes[:meanChanges600gens]); digits=7))
print("\nmean env Robustness before evolution when 1224 env changes = ", round(mean(resultsDataFrame1224changes[:meanChangesBeforeEvol]); digits=7))
print("\nmean env Robustness after 100 gens when 1224 env changes = ", round(mean(resultsDataFrame1224changes[:meanChanges100gens]); digits=7))
print("\nmean env Robustness after 600 gens when 1224 env changes = ", round(mean(resultsDataFrame1224changes[:meanChanges600gens]); digits=7))
# Print out stds of the trials for each case:
print("\n\nstd env Robustness before evolution when 100 env changes = ", round(std(resultsDataFrame100changes[:meanChangesBeforeEvol]); digits=7))
print("\nstd env Robustness after 100 gens when 100 env changes = ", round(std(resultsDataFrame100changes[:meanChanges100gens]); digits=7))
print("\nstd env Robustness after 600 gens when 100 env changes = ", round(std(resultsDataFrame100changes[:meanChanges600gens]); digits=7))
print("\nstd env Robustness before evolution when 1224 env changes = ", round(std(resultsDataFrame1224changes[:meanChangesBeforeEvol]); digits=7))
print("\nstd env Robustness after 100 gens when 1224 env changes = ", round(std(resultsDataFrame1224changes[:meanChanges100gens]); digits=7))
print("\nstd env Robustness after 600 gens when 1224 env changes = ", round(std(resultsDataFrame1224changes[:meanChanges600gens]); digits=7))

# Plot comparison of 100 vs 1224 envchanges including number of stable perturbations:
## Before evolution:
plot(stablePerturbationsResults100changes[:stablePerturbationsBeforeEvol]./100,resultsDataFrame100changes[:meanChangesBeforeEvol],seriestype=:scatter,title="100 envChanges, Before Evol", label="100",ylabel="env Robustness")
plot100vs1224StablePerturbationsBeforeEvol = plot!(stablePerturbationsResults1224changes[:stablePerturbationsBeforeEvol]./1224,resultsDataFrame1224changes[:meanChangesBeforeEvol],seriestype=:scatter, title="100 vs 1224 envChanges, Before Evol", label="1224",ylabel="env Robustness")
## After 100 gens:
plot(stablePerturbationsResults100changes[:stablePerturbations100gens]./100,resultsDataFrame100changes[:meanChanges100gens],seriestype=:scatter, title="100 envChanges, 100 generations", label="100")
plot100vs1224StablePerturbations100gens = plot!(stablePerturbationsResults1224changes[:stablePerturbations100gens]./1224,resultsDataFrame1224changes[:meanChanges100gens],seriestype=:scatter, title="100 vs 1224 envChanges, 100 generations", label="1224",ylabel="env Robustness")
## After 600 gens:
plot(stablePerturbationsResults100changes[:stablePerturbations600gens]./100,resultsDataFrame100changes[:meanChanges600gens],seriestype=:scatter, title="100 envChanges, 600 generations", label="100")
plot100vs1224StablePerturbations600gens = plot!(stablePerturbationsResults1224changes[:stablePerturbations600gens]./1224,resultsDataFrame1224changes[:meanChanges600gens],seriestype=:scatter, title="100 vs 1224 envChanges, 600 generations", label="1224",xlabel="% stable perturbations",ylabel="env Robustness")
plot(plot100vs1224StablePerturbationsBeforeEvol,plot100vs1224StablePerturbations100gens,plot100vs1224StablePerturbations600gens,layout=(3,1))
print("\nmean stable perturbations before evolution when 100 env changes = ", round(mean(stablePerturbationsResults100changes[:stablePerturbationsBeforeEvol])/100; digits=7))
print("\nmean stable perturbations after 100 gens when 100 env changes = ", round(mean(stablePerturbationsResults100changes[:stablePerturbations100gens])/100; digits=7))
print("\nmean stable perturbations after 600 gens when 100 env changes = ", round(mean(stablePerturbationsResults100changes[:stablePerturbations600gens])/100; digits=7))
print("\nmean stable perturbations before evolution when 1224 env changes = ", round(mean(stablePerturbationsResults1224changes[:stablePerturbationsBeforeEvol])/1224; digits=7))
print("\nmean stable perturbations after 100 gens when 1224 env changes = ", round(mean(stablePerturbationsResults1224changes[:stablePerturbations100gens])/1224; digits=7))
print("\nmean stable perturbations after 600 gens when 1224 env changes = ", round(mean(stablePerturbationsResults1224changes[:stablePerturbations600gens])/1224; digits=7))

# Plot comparison of envchanges throughout evolution:
means100changes = [mean(resultsDataFrame100changes[:meanChangesBeforeEvol]),mean(resultsDataFrame100changes[:meanChanges100gens]),mean(resultsDataFrame100changes[:meanChanges600gens])]
means1224changes = [mean(resultsDataFrame1224changes[:meanChangesBeforeEvol]),mean(resultsDataFrame1224changes[:meanChanges100gens]),mean(resultsDataFrame1224changes[:meanChanges600gens])]
plot([[0,100,600],[0,100,600]],[means100changes,means1224changes],title="env robustness throughout evolution for multiple different trials",label=["100 changes","1224 changes"],ylabel="env robustness",xlabel="generations")
# --> means of different populations are very close in trend, so I don't think 100 vs 1224 envchanges matters for our analysis of multiple different populations when measuring environmental robustness behavior

## Test with polycomb
POLYCOMBFLAG = true
include("individuals.jl")
# found env. robustness increases, had lab meeting and decided env robustness might not be good measure to focus on because we do not know the size of the env. perturbation like how far pushing away from original environment.

## Test with no selection but must compare to individual before perturbation instead of founder
###############################################################################



################################################################################
# Test when WxE set to zero and SAURABVERSION is true, they should be the same.
# Start with the same founder also.
################################################################################
## Set environment interaction matrix to zero and see if get env robustness trend
## like see for Saurabh version:
# Initialize number of trials to run for different populations to compare:
numTests = 1
# Initialize dataframes to store results:
resultsDataFrame100changes = DataFrame(meanChangesBeforeEvol = Float64[], meanChanges100gens = Float64[], meanChanges600gens = Float64[])
stdResultsDataFrame100changes = DataFrame(stdChangesBeforeEvol = Float64[], stdChanges100gens = Float64[], stdChanges600gens = Float64[])
stablePerturbationsResults100changes = DataFrame(stablePerturbationsBeforeEvol = Float64[], stablePerturbations100gens = Float64[], stablePerturbations600gens = Float64[])
resultsDataFrame100changesSaurabhV = DataFrame(meanChangesBeforeEvol = Float64[], meanChanges100gens = Float64[], meanChanges600gens = Float64[])
stdResultsDataFrame100changesSaurabhV = DataFrame(stdChangesBeforeEvol = Float64[], stdChanges100gens = Float64[], stdChanges600gens = Float64[])
stablePerturbationsResults100changesSaurabhV = DataFrame(stablePerturbationsBeforeEvol = Float64[], stablePerturbations100gens = Float64[], stablePerturbations600gens = Float64[])
# Run test:
for i = 1:numTests
    # Do with Saurabh version of population:
    print("\nmeasure ", i, " population for Saurabh version")
    SAURABVERSION = true
    include("types.jl")
    include("individuals.jl")
    include("population.jl")
    saurabhpop = genpop()
    # before evolution:
    founderEnvRobust(saurabhpop.founder)
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(saurabhpop.founder)
    stablePerturbationsChangesBeforeEvol = map(x->environmentalRobustness(saurabhpop.individuals[x], mutatedEnvsMatrix, saurabhpop.founder.envRobustness), 1:N)
    individualEnvRobustnessChangesBeforeEvol = map(x->saurabhpop.individuals[x].envRobustness, 1:N)
    # after 100 generations of evolution:
    for j = 1:100
        update(saurabhpop)
    end
    print("\nmeasure after 100 gens:")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(saurabhpop.founder)
    stablePerturbationsChanges100gens = map(x->environmentalRobustness(saurabhpop.individuals[x], mutatedEnvsMatrix, saurabhpop.founder.envRobustness), 1:N)
    individualEnvRobustnessChanges100gens = map(x->saurabhpop.individuals[x].envRobustness, 1:N)
    # after 600 generations of evolution:
    for j = 1:500
        update(saurabhpop)
    end
    print("\nmeasure after 600 gens:")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(saurabhpop.founder)
    stablePerturbationsChanges600gens = map(x->environmentalRobustness(saurabhpop.individuals[x], mutatedEnvsMatrix, saurabhpop.founder.envRobustness), 1:N)
    individualEnvRobustnessChanges600gens = map(x->saurabhpop.individuals[x].envRobustness, 1:N)
    # Add all the results to dataframes:
    push!(resultsDataFrame100changesSaurabhV, [mean(individualEnvRobustnessChangesBeforeEvol),mean(individualEnvRobustnessChanges100gens),mean(individualEnvRobustnessChanges600gens)])
    push!(stdResultsDataFrame100changesSaurabhV, [std(individualEnvRobustnessChangesBeforeEvol),std(individualEnvRobustnessChanges100gens),std(individualEnvRobustnessChanges600gens)])
    push!(stablePerturbationsResults100changesSaurabhV, [mean(stablePerturbationsChangesBeforeEvol),mean(stablePerturbationsChanges100gens),mean(stablePerturbationsChanges600gens)])

    # Now with environment
    print("\nmeasure ", i, " population")
    const WMAT = saurabhpop.founder.network # will always work because this individual is already stable; if give random W matrix as WMAT then genpop() may not be able to get a stable founder --> error will return if so
    const INP = saurabhpop.founder.initstate
    const OPT = (Float64)[]
    SAURABVERSION = false
    include("types.jl")
    include("individuals.jl")
    include("population.jl")
    pop = genpop()
    pop.founder.EnvIndInter = zeros(50,50)
    map(x->pop.individuals[x].EnvIndInter = zeros(50,50),(1:N));
    # before evolution:
    founderEnvRobust(pop.founder)
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbationsChangesBeforeEvol = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustnessChangesBeforeEvol = map(x->pop.individuals[x].envRobustness, 1:N)
    # after 100 generations of evolution:
    for j = 1:100
        update(pop)
    end
    print("\nmeasure after 100 gens:")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbationsChanges100gens = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustnessChanges100gens = map(x->pop.individuals[x].envRobustness, 1:N)
    # after 600 generations of evolution:
    for j = 1:500
        update(pop)
    end
    print("\nmeasure after 500 gens:")
    mutatedEnvsMatrix = generateMutatedEnvsMatrix(pop.founder)
    stablePerturbationsChanges600gens = map(x->environmentalRobustness(pop.individuals[x], mutatedEnvsMatrix, pop.founder.envRobustness), 1:N)
    individualEnvRobustnessChanges600gens = map(x->pop.individuals[x].envRobustness, 1:N)
    # Add all the results to dataframes:
    push!(resultsDataFrame100changes, [mean(individualEnvRobustnessChangesBeforeEvol),mean(individualEnvRobustnessChanges100gens),mean(individualEnvRobustnessChanges600gens)])
    push!(stdResultsDataFrame100changes, [std(individualEnvRobustnessChangesBeforeEvol),std(individualEnvRobustnessChanges100gens),std(individualEnvRobustnessChanges600gens)])
    push!(stablePerturbationsResults100changes, [mean(stablePerturbationsChangesBeforeEvol),mean(stablePerturbationsChanges100gens),mean(stablePerturbationsChanges600gens)])
end
# Plot results:
meansEnvInteractionZero = [mean(resultsDataFrame100changes[:meanChangesBeforeEvol]),mean(resultsDataFrame100changes[:meanChanges100gens]),mean(resultsDataFrame100changes[:meanChanges600gens])]
# meansSaurabh = [mean(resultsDataFrame100changesSaurabhV[:meanChangesBeforeEvol]),mean(resultsDataFrame100changesSaurabhV[:meanChanges100gens]),mean(resultsDataFrame100changesSaurabhV[:meanChanges600gens])]
# plot([[0,100,600],[0,100,600]], [meansEnvInteractionZero,meansSaurabh], title="env robustness during evolution for  env inter 0 vs. saurabh v", label=["EnvInter0","SaurabhV"],ylabel="env robustness",xlabel="generations")
plot([0,100,600],meansEnvInteractionZero,title="mean env robustness when env interaction matrix 0")
################################################################################
