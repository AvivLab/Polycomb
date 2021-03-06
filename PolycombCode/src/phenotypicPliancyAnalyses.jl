# Created by Maryl Lambros on 09/11/2019
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
# indir = joinpath("..","input")
# indir = "input"
indir = joinpath("/polycomb", "input")
constantsFile = joinpath(indir, configfile) # constants.jl used is in the input folder in julia folder
include(constantsFile)

### Added by Aviv
#if isempty(ARGS);
#else
#	MUTRATE=parse(Float64, ARGS[1])
#	POLYMUTRATE=parse(Float64, ARGS[1])
#	ENVMUTRATE=parse(Float64, ARGS[1])
#	C=parse(Float64, ARGS[2])
#	SELSTR=parse(Float64, ARGS[3])
#	SIGSTR=parse(Float64, ARGS[4])
#	GENETHRESH=parse(Float64, ARGS[5])
#end
### End Added by Aviv

print("Running script...")

include("utilities.jl")
include("types.jl")
include("individuals.jl")
include("population.jl")
include("measure.jl")
include("ipaNetworkFileConversion.jl")

#-------------------------------------------------------------------------
# Functions to run phenotypic pliancy analyses
#-------------------------------------------------------------------------
function meanDistanceBtwPops(pop1ToTest::Population, pop2ToTest::Population, envState1ToTest::Int64, envState2ToTest::Int64)
    distanceBtwPops = map(x -> if (pop1ToTest.individuals[x].stable == trues(2)) & (pop2ToTest.individuals[x].stable == trues(2)); return sum(abs.(pop1ToTest.individuals[x].develstate[:,envState1ToTest] - pop2ToTest.individuals[x].develstate[:,envState2ToTest]))/G; end, 1:N)
    distanceBtwPops = convert(Array{Float64}, distanceBtwPops[distanceBtwPops.!=nothing])
    return mean(distanceBtwPops)
end

function meanDistanceToOptStateBtwPops(pop1ToTest::Population, pop2ToTest::Population, envState1ToTest::Int64, envState2ToTest::Int64)
    distanceBtwPops = map(x -> if (pop1ToTest.individuals[x].stable == trues(2)) & (pop2ToTest.individuals[x].stable == trues(2)); return sum(abs.(pop1ToTest.individuals[x].develstate[:,envState1ToTest] - pop2ToTest.individuals[x].optstate[:,envState2ToTest]))/G; end, 1:N)
    distanceBtwPops = convert(Array{Float64}, distanceBtwPops[distanceBtwPops.!=nothing])
    return mean(distanceBtwPops)
end

function phenotypicPliancyPercentage(pop1ToTest::Population, pop2ToTest::Population, envState1ToTest::Int64, envState2ToTest::Int64)
    distanceBtwPops = map(x -> if (pop1ToTest.individuals[x].stable == trues(2)) & (pop2ToTest.individuals[x].stable == trues(2)); return sum(abs.(pop1ToTest.individuals[x].develstate[:,envState1ToTest] - pop2ToTest.individuals[x].develstate[:,envState2ToTest]))/G; end, 1:N)
    distanceBtwPops = convert(Array{Float64}, distanceBtwPops[distanceBtwPops.!=nothing])
    return round(length(findall(y -> y > MINDISTFORPLIANCY, distanceBtwPops))/length(distanceBtwPops)*100, digits = 2)
end

function meanOfResults(distanceResults::Array{Float64})
    distanceResultsNaNsRemoved = removeNaNs(distanceResults)
    return round(mean(distanceResultsNaNsRemoved), digits = 4)
end

function removeNaNs(distanceResults::Array{Float64})
    return distanceResults[findall(!isnan,distanceResults)]
end

function meanRatioPcgBreakageToDistance(popWithBroken::Population, popWithIntact::Population, envStateToTest::Int64)
    distanceBtwPops = map(x -> if (popWithBroken.individuals[x].stable == trues(2)) & (popWithIntact.individuals[x].stable == trues(2)); return (length(findall(y->y == 0, popWithIntact.individuals[x].polycombstate[:,envStateToTest]))/ENVS)/(sum(abs.(popWithBroken.individuals[x].develstate[:,envStateToTest] - popWithIntact.individuals[x].develstate[:,envStateToTest]))/G); end, 1:N)
    distanceBtwPops = convert(Array{Float64}, distanceBtwPops[distanceBtwPops.!=nothing])
    return mean(distanceBtwPops)
end

function gatherDevelStateResultsRightAfterEvol(pop::Population,envstate1::Int64,envstate2::Int64)
    y1 = zeros(Float64, G, N);
    y2 = zeros(Float64, G, N);
    for a = 1:N
            y1[:,a] = pop.individuals[a].develstate[:,envstate1];
            y2[:,a] = pop.individuals[a].develstate[:,envstate2];
    end
    return y1, y2
end

function gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcs::Population,popWithIntactPrcs::Population,envstate::Int64)
    y1 = zeros(Float64, G, N);
    y2 = zeros(Float64, G, N);
    for a = 1:N
            y1[:,a] = popWithBrokenPrcs.individuals[a].develstate[:,envstate];
            y2[:,a] = popWithIntactPrcs.individuals[a].develstate[:,envstate];
    end
    return y1, y2
end

function gatherPolycombStateResultsRightAfterEvol(pop::Population,envstate1::Int64,envstate2::Int64)
    y1 = zeros(Float64, G, N);
    y2 = zeros(Float64, G, N);
    for a = 1:N
            y1[:,a] = pop.individuals[a].polycombstate[:,envstate1];
            y2[:,a] = pop.individuals[a].polycombstate[:,envstate2];
    end
    return y1, y2
end

function gatherPolycombVecResultsRightAfterEvol(pop::Population,envstate1::Int64,envstate2::Int64)
    y1 = zeros(Float64, G, N);
    y2 = zeros(Float64, G, N);
    for a = 1:N
            y1[:,a] = pop.individuals[a].polycombvec[:,envstate1];
            y2[:,a] = pop.individuals[a].polycombvec[:,envstate2];
    end
    return y1, y2
end

function savingGeneExpResults(geneExpResults::Array{Float64,2}, namesForDf::Array{Symbol,1}, filenameForSaving::String, directoryToSaveIn::String)
    geneExpResultsDf = convert(DataFrame, geneExpResults);
    names!(geneExpResultsDf, namesForDf);
    CSV.write(joinpath(directoryToSaveIn, filenameForSaving), geneExpResultsDf)
end


##########################################################################
#-------------------------------------------------------------------------
# Run phenotypic pliancy analyses
#-------------------------------------------------------------------------
#include(constantsFile)
if isempty(ARGS)
	timestamp = gentimestamp()
else
	timestamp=ARGS[6]
end

# dataIntactAndBrokenoutdir = joinpath("..", "results")
dataIntactAndBrokenoutdir = "results"
dataIntactAndBrokensimdir = joinpath(dataIntactAndBrokenoutdir, timestamp)
# save constants.jl file used for run of this code
configpath = joinpath(dataIntactAndBrokensimdir, configfile)
if !isdir(dataIntactAndBrokensimdir)
    mkdir(dataIntactAndBrokensimdir)
end
cp(constantsFile, configpath)

if isempty(ARGS);
else
    MUTRATE=parse(Float64, ARGS[1])
    POLYMUTRATE=parse(Float64, ARGS[1])
    ENVMUTRATE=parse(Float64, ARGS[1])
    C=parse(Float64, ARGS[2])
    SELSTR=parse(Float64, ARGS[3])
    SIGSTR=parse(Float64, ARGS[4])
    GENETHRESH=parse(Float64, ARGS[5])

	f=open(configpath,"a")
	println(f,"\nMUTRATE = ", MUTRATE)
	println(f,"POLYMUTRATE = ", POLYMUTRATE)
	println(f,"ENVMUTRATE = ", ENVMUTRATE)
	println(f,"C = ", C)
	println(f,"SELSTRE = ", SELSTR)
	println(f,"SIGSTR = ", SIGSTR)
	println(f,"GENETHRESH = ", GENETHRESH)
	#println(f,"timestamp = ", timestamp, "\n")
	close(f)
end
print("\n",dataIntactAndBrokensimdir)
for h = 1:1 # when 1:1 then break BOTH PRCs only; if 3 then break each individually and then combination of two
    prcsToRemove = [[1,2]][h];
    #prcsToRemove = [[1],[2],[1,2]][h];
    if length(prcsToRemove) == 1 # if only removing one PRC
        prcsRemovedFolderName = string("PrcRemoved",prcsToRemove[1]);
    elseif length(prcsToRemove) > 1 # if removing more than one PRC
        prcsRemovedFolderName = string("PrcRemoved",join(map(x->string(prcsToRemove[x]),1:length(prcsToRemove)),"and"));
    else # if no PRC to remove
        error(string("Error: No PRC selected to remove."))
    end
    folderForDataStorage = joinpath(dataIntactAndBrokensimdir, prcsRemovedFolderName)
    if !isdir(folderForDataStorage)
        mkdir(folderForDataStorage)
    end
    ## Variables to store results when use initstate as starting point of iterations when testing phenotype pliancy:
    # distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenStayInEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv1WithIntactPrcWhenStayInEnv1 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenStayInEnv2 = zeros(Float64,INDTRIALS);
    # distanceToOriginalDevelStateInGivenEnv2WithIntactPrcWhenStayInEnv2 = zeros(Float64,INDTRIALS);
    # distanceWithNewEnvResultsReplaceEnv1 = zeros(Float64,INDTRIALS);
    # distanceWithNewEnvResultsReplaceEnv2 = zeros(Float64,INDTRIALS);
    # ratioOfPcgBreakageToDistanceMoveInNewEnvResultsWhenReplaceEnv1 = zeros(Float64,INDTRIALS);
    # ratioOfPcgBreakageToDistanceMoveInNewEnvResultsWhenReplaceEnv2 = zeros(Float64,INDTRIALS);
    # distanceWithNewEnvForEnv1WhenPrcIntactResults = zeros(Float64,INDTRIALS);
    # distanceWithNewEnvForEnv2WhenPrcIntactResults = zeros(Float64,INDTRIALS);
    # distanceWithNewEnvForEnv1WhenPrcBrokenResults = zeros(Float64,INDTRIALS);
    # distanceWithNewEnvForEnv2WhenPrcBrokenResults = zeros(Float64,INDTRIALS);
    ## Variables to store results
    numStableIndWhenPrcIntactAndNewEnvReplaceEnv1 = zeros(Float64,INDTRIALS);
    numStableIndWhenPrcBrokenAndNewEnvReplaceEnv1 = zeros(Float64,INDTRIALS);
    numStableIndWhenPrcIntactAndNewEnvReplaceEnv2 = zeros(Float64,INDTRIALS);
    numStableIndWhenPrcBrokenAndNewEnvReplaceEnv2 = zeros(Float64,INDTRIALS);
    numStableIndWhen1PrcElementBrokenAndNewEnvReplaceEnv1 = zeros(Float64,INDTRIALS);
    numStableIndWhen1PrcElementBrokenAndNewEnvReplaceEnv2 = zeros(Float64,INDTRIALS);
    fitnessOfIndividualsAtStartOfEvolution = zeros(Float64, 2, INDTRIALS);
    fitnessOfIndividualsAtEndOfEvolution = zeros(Float64, 2, INDTRIALS);
    overallFitnessBeforeEvolution = zeros(Float64,INDTRIALS);
    overallFitnessAfterEvolution = zeros(Float64,INDTRIALS);
    averageNumOfGenesSuppressedByPrcs = zeros(Float64, 2, INDTRIALS);

    ## Variables to store results when use final develstate as starting point of iterations when testing for phenotypic pliancy instead of initstate:
    distanceWhenPolycombIntactInEnv1 = zeros(Float64,INDTRIALS);
    distanceWhenPolycombIntactInEnv2 = zeros(Float64,INDTRIALS);
    distanceWhenPolycombBrokenInEnv1 = zeros(Float64,INDTRIALS);
    distanceWhenPolycombBrokenInEnv2 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv2 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv2 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv1 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv1 = zeros(Float64,INDTRIALS);
    distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2 = zeros(Float64,INDTRIALS);
    distanceWhenPolycombIntactInEnv1AndNewEnv = zeros(Float64,INDTRIALS);
    distanceWhenPolycombBrokenInEnv1AndNewEnv = zeros(Float64,INDTRIALS);
    distanceWhenPolycombIntactInEnv2AndNewEnv = zeros(Float64,INDTRIALS);
    distanceWhenPolycombBrokenInEnv2AndNewEnv = zeros(Float64,INDTRIALS);

    ## Variables for measuring percent phenotypically pliant
    percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1 = zeros(Float64,INDTRIALS);
    # percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv2 = zeros(Float64,INDTRIALS);
    percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1 = zeros(Float64,INDTRIALS);
    # percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv2 = zeros(Float64,INDTRIALS);
    # percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv1 = zeros(Float64,INDTRIALS);
    percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2 = zeros(Float64,INDTRIALS);
    # percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv1 = zeros(Float64,INDTRIALS);
    percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2 = zeros(Float64,INDTRIALS);
    percentPliantWhenPolycombIntactInEnv1AndNewEnv = zeros(Float64,INDTRIALS);
    percentPliantWhenPolycombBrokenInEnv1AndNewEnv = zeros(Float64,INDTRIALS);
    percentPliantWhenPolycombIntactInEnv2AndNewEnv = zeros(Float64,INDTRIALS);
    percentPliantWhenPolycombBrokenInEnv2AndNewEnv = zeros(Float64,INDTRIALS);

    newEnv3StatesToPrint = zeros(Float64, G, INDTRIALS);

    polycombStateAfterEvolEnv1 = zeros(Float64, G, N*INDTRIALS);
    polycombStateAfterEvolEnv2 = zeros(Float64, G, N*INDTRIALS);
    polycombVecAfterEvolEnv1 = zeros(Float64, G, N*INDTRIALS);
    polycombVecAfterEvolEnv2 = zeros(Float64, G, N*INDTRIALS);

    geneExpAfterEvolEnv1 = zeros(Float64, G, N*INDTRIALS);
    geneExpAfterEvolEnv2 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgIntactEnv1 = zeros(Float64, G, N*INDTRIALS); # For grant save results when start with initstate
    geneExpResultsPcgIntactEnv2 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgBrokenEnv1 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgBrokenEnv2 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgIntactEnv1ReplaceEnv2 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgIntactEnv2ReplaceEnv1 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgBrokenEnv1ReplaceEnv2 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsPcgBrokenEnv2ReplaceEnv1 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv1 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv1 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv2 = zeros(Float64, G, N*INDTRIALS);
    geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv2 = zeros(Float64, G, N*INDTRIALS);
    if SAMESTARTINGPOPFORINDTRIALSFLAG == true
        Random.seed!(RANDSEEDNUM[1]);
        pop = genpop();
    end
    for j = 1:INDTRIALS
        print("\nevolution for pop. number: ", j, " gens = ", GENS, " and PRC to remove = ", prcsToRemove)
        Random.seed!(RANDSEEDNUM[j+1]);
        if SAMESTARTINGPOPFORINDTRIALSFLAG == false
            pop1 = genpop();
        else
            pop1 = deepcopy(pop);
        end
        fitnessOfIndividualsAtStartOfEvolution[:,j] = mean(map(x -> pop1.individuals[x].fitnessUnderEachEnv, 1:N))
        overallFitnessBeforeEvolution[j] = mean(map(x -> pop1.individuals[x].fitness, 1:N))
        for i = 1:GENS
            update(pop1);
        end
        geneExpAfterEvolEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpAfterEvolEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsRightAfterEvol(pop1, 1, 2);
        fitnessOfIndividualsAtEndOfEvolution[:,j] = mean(map(x -> pop1.individuals[x].fitnessUnderEachEnv, 1:N));
        overallFitnessAfterEvolution[j] = mean(map(x -> pop1.individuals[x].fitness, 1:N));
        # Find polycomb mechanism results after evolution:
        averageNumOfGenesSuppressedByPrcs[:,j] = mean(map(z -> map(y -> length(findall(x -> x == 0, pop1.individuals[z].polycombstate[:,y])), 1:2), 1:N))
        # Find to save polycomb mechanism suppression information for all genes:
        polycombStateAfterEvolEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], polycombStateAfterEvolEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherPolycombStateResultsRightAfterEvol(pop1,1,2);
        polycombVecAfterEvolEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], polycombVecAfterEvolEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherPolycombVecResultsRightAfterEvol(pop1,1,2);
        popWithBrokenPrcs = deepcopy(pop1);
        popWithoutBrokenPrcs = deepcopy(pop1);
        if BREAKENTIREPRC != true # If break only 1 target element when only 1 PRC:
            map(x -> breakpoly(popWithBrokenPrcs.individuals[x], 2), 1:N)
        else # If multiple PRCs:
            map(x -> breakPRC(popWithBrokenPrcs.individuals[x], prcsToRemove, 1), 1:N) # break PRC (prcsToRemove) in env 1
            map(x -> breakPRC(popWithBrokenPrcs.individuals[x], prcsToRemove, 2), 1:N) # break PRC (prcsToRemove) in env 2
        end
        # Populations to test in totally different environment than ones evolved in:
        popWithBrokenPrcsAndTestNewEnv = deepcopy(popWithBrokenPrcs);
        popWithoutBrokenPrcsAndTestNewEnv = deepcopy(popWithoutBrokenPrcs);
        # Populations to test for stability when start with develstate instead of initstate:
        popWithBrokenPolycomb = deepcopy(popWithBrokenPrcs);
        popWithIntactPolycomb = deepcopy(popWithoutBrokenPrcs);
        # Populations to test stability when start with develstate instead of initstate and replace env2 with env1 and vice versus:
        popWithBrokenPolycombReplaceEnvs = deepcopy(popWithBrokenPrcs);
        popWithIntactPolycombReplaceEnvs = deepcopy(popWithoutBrokenPrcs);

        ### Test for stability when start with develstate instead of initstate:
        # In Env 1:
        map(x -> iterateindForStability(popWithIntactPolycomb.individuals[x], 1, pop1.founder.EnvState[:,1], 1, pop1.individuals[x]), 1:N)
        map(x -> iterateindForStability(popWithBrokenPolycomb.individuals[x], 1, pop1.founder.EnvState[:,1], 1, pop1.individuals[x]), 1:N)
        distanceWhenPolycombIntactInEnv1[j] = meanDistanceBtwPops(popWithIntactPolycomb, pop1, 1, 1);
        distanceWhenPolycombBrokenInEnv1[j] = meanDistanceBtwPops(popWithBrokenPolycomb, pop1, 1, 1);
        # Save results:
        geneExpResultsPcgBrokenEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPolycomb,popWithIntactPolycomb,1)
        # In Env 2:
        map(x -> iterateindForStability(popWithIntactPolycomb.individuals[x], 2, pop1.founder.EnvState[:,2], 2, pop1.individuals[x]), 1:N)
        map(x -> iterateindForStability(popWithBrokenPolycomb.individuals[x], 2, pop1.founder.EnvState[:,2], 2, pop1.individuals[x]), 1:N)
        distanceWhenPolycombIntactInEnv2[j] = meanDistanceBtwPops(popWithIntactPolycomb, pop1, 2, 2);
        distanceWhenPolycombBrokenInEnv2[j] = meanDistanceBtwPops(popWithBrokenPolycomb, pop1, 2, 2);
        # Save results:
        geneExpResultsPcgBrokenEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPolycomb,popWithIntactPolycomb,2)
        ## If evolve in env 1 (develstate from env 1) but place in env 2:
        map(x -> iterateindForStability(popWithIntactPolycombReplaceEnvs.individuals[x], 2, pop1.founder.EnvState[:,2], 1, pop1.individuals[x]), 1:N)
        map(x -> iterateindForStability(popWithBrokenPolycombReplaceEnvs.individuals[x], 2, pop1.founder.EnvState[:,2], 1, pop1.individuals[x]), 1:N)
        # compared to env 1:
        distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1[j] = meanDistanceBtwPops(popWithIntactPolycombReplaceEnvs, pop1, 2, 1);
        pop1ToTest= deepcopy(popWithIntactPolycombReplaceEnvs)
        pop2ToTest=deepcopy(pop1)
        envState1ToTest=2
        envState2ToTest=2
        distanceBtwPops = map(x -> if (pop1ToTest.individuals[x].stable == trues(2)) & (pop2ToTest.individuals[x].stable == trues(2)); return sum(abs.(pop1ToTest.individuals[x].develstate[:,envState1ToTest] - pop2ToTest.individuals[x].develstate[:,envState2ToTest]))/G; end, 1:N)
        distanceBtwPops = convert(Array{Float64}, distanceBtwPops[distanceBtwPops.!=nothing])
        distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1[j] = meanDistanceBtwPops(popWithBrokenPolycombReplaceEnvs, pop1, 2, 1);
        percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1[j] = phenotypicPliancyPercentage(popWithIntactPolycombReplaceEnvs, pop1, 2, 1);
        percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1[j] = phenotypicPliancyPercentage(popWithBrokenPolycombReplaceEnvs, pop1, 2, 1);
        # compared to env 2:
        distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv2[j] = meanDistanceBtwPops(popWithIntactPolycombReplaceEnvs, pop1, 2, 2);
        distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv2[j] = meanDistanceBtwPops(popWithBrokenPolycombReplaceEnvs, pop1, 2, 2);
        # percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv2[j] = phenotypicPliancyPercentage(popWithIntactPolycombReplaceEnvs, pop1, 2, 2);
        # percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv2[j] = phenotypicPliancyPercentage(popWithBrokenPolycombReplaceEnvs, pop1, 2, 2);
        # Save results:
        geneExpResultsPcgBrokenEnv2ReplaceEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv2ReplaceEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPolycombReplaceEnvs,popWithIntactPolycombReplaceEnvs,2)
        ## If evolve in env 2 (develstate from env 1) but place in env 1:
        map(x -> iterateindForStability(popWithIntactPolycombReplaceEnvs.individuals[x], 1, pop1.founder.EnvState[:,1], 2, pop1.individuals[x]), 1:N)
        map(x -> iterateindForStability(popWithBrokenPolycombReplaceEnvs.individuals[x], 1, pop1.founder.EnvState[:,1], 2, pop1.individuals[x]), 1:N)
        # compared to env 1:
        distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv1[j] = meanDistanceBtwPops(popWithIntactPolycombReplaceEnvs, pop1, 1, 1);
        distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv1[j] = meanDistanceBtwPops(popWithBrokenPolycombReplaceEnvs, pop1, 1, 1);
        # percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv1[j] = phenotypicPliancyPercentage(popWithIntactPolycombReplaceEnvs, pop1, 1, 1);
        # percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv1[j] = phenotypicPliancyPercentage(popWithBrokenPolycombReplaceEnvs, pop1, 1, 1);
        # compared to env 2:
        distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2[j] = meanDistanceBtwPops(popWithIntactPolycombReplaceEnvs, pop1, 1, 2);
        distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2[j] = meanDistanceBtwPops(popWithBrokenPolycombReplaceEnvs, pop1, 1, 2);
        percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2[j] = phenotypicPliancyPercentage(popWithIntactPolycombReplaceEnvs, pop1, 1, 2);
        percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2[j] = phenotypicPliancyPercentage(popWithBrokenPolycombReplaceEnvs, pop1, 1, 2);
        # Save results:
        geneExpResultsPcgBrokenEnv1ReplaceEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv1ReplaceEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPolycombReplaceEnvs,popWithIntactPolycombReplaceEnvs,1)
        ### Start at initstate instead of develstate:
        ## Distance when stay in environment 1 and 2, respectively, after breaking or keeping PRC intact:
        # Env 1:
        # map(x -> iterateind2(popWithBrokenPrcs.individuals[x], 1, popWithoutBrokenPrcs.founder.EnvState[:,1]), 1:N) # test if stay in 1st environmental state when break PRC
        # map(x -> iterateind2(popWithoutBrokenPrcs.individuals[x], 1, popWithoutBrokenPrcs.founder.EnvState[:,1]), 1:N) # test if stay in 1st environmental state when PRCs intact
        # # Save results:
        # geneExpResultsPcgBrokenEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcs,popWithoutBrokenPrcs,1)
        # # Measurements:
        # distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenStayInEnv1[j] = meanDistanceBtwPops(popWithBrokenPrcs, pop1, 1, 1)
        # distanceToOriginalDevelStateInGivenEnv1WithIntactPrcWhenStayInEnv1[j] = meanDistanceBtwPops(popWithoutBrokenPrcs, pop1, 1, 1)
        # # In env 2:
        # map(x -> iterateind2(popWithBrokenPrcs.individuals[x], 2, popWithoutBrokenPrcs.founder.EnvState[:,2]), 1:N) # test if stay in 1st environmental state when break PRC
        # map(x -> iterateind2(popWithoutBrokenPrcs.individuals[x], 2, popWithoutBrokenPrcs.founder.EnvState[:,2]), 1:N) # test if stay in 1st environmental state when PRCs intact
        # # Save results:
        # geneExpResultsPcgBrokenEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcs,popWithoutBrokenPrcs,2)
        # # Measurements:
        # distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenStayInEnv2[j] = meanDistanceBtwPops(popWithBrokenPrcs, pop1, 2, 2)
        # distanceToOriginalDevelStateInGivenEnv2WithIntactPrcWhenStayInEnv2[j] = meanDistanceBtwPops(popWithoutBrokenPrcs, pop1, 2, 2)
        # ## Distance with broken PRC when replace env 1 with env 2:
        # map(x -> iterateind2(popWithBrokenPrcs.individuals[x], 1, popWithoutBrokenPrcs.founder.EnvState[:,2]), 1:N) # test if change to 1st environmental state when break PRC
        # map(x -> iterateind2(popWithoutBrokenPrcs.individuals[x], 1, popWithoutBrokenPrcs.founder.EnvState[:,2]), 1:N) # test if change to 1st environmental state when PRCs intact
        # # Save results:
        # geneExpResultsPcgBrokenEnv2ReplaceEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv2ReplaceEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcs,popWithoutBrokenPrcs,1)
        # # Measurements:
        # distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceBtwPops(popWithBrokenPrcs, pop1, 1, 2)
        # distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceBtwPops(popWithBrokenPrcs, pop1, 1, 1)
        # distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceToOptStateBtwPops(popWithBrokenPrcs, pop1, 1, 2)
        # distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceToOptStateBtwPops(popWithBrokenPrcs, pop1, 1, 1)
        # # Distance without broken PRC:
        # distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceBtwPops(popWithoutBrokenPrcs, pop1, 1, 2)
        # distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceBtwPops(popWithoutBrokenPrcs, pop1, 1, 1)
        # distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceToOptStateBtwPops(popWithoutBrokenPrcs, pop1, 1, 2)
        # distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2[j] = meanDistanceToOptStateBtwPops(popWithoutBrokenPrcs, pop1, 1, 1)
        # ## Distance with broken PRC when replace env 2 with env 1:
        # # Distance with broken PRC:
        # map(x -> iterateind2(popWithBrokenPrcs.individuals[x], 2, popWithoutBrokenPrcs.founder.EnvState[:,1]), 1:N) # test if change to 2nd environmental state when break PRC
        # map(x -> iterateind2(popWithoutBrokenPrcs.individuals[x], 2, popWithoutBrokenPrcs.founder.EnvState[:,1]), 1:N) # test if change to 2nd environmental state when PRCs intact
        # # Save results:
        # geneExpResultsPcgBrokenEnv1ReplaceEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsPcgIntactEnv1ReplaceEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcs,popWithoutBrokenPrcs,2)
        # distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceBtwPops(popWithBrokenPrcs, pop1, 2, 1)
        # distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceBtwPops(popWithBrokenPrcs, pop1, 2, 2)
        # distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceToOptStateBtwPops(popWithBrokenPrcs, pop1, 2, 1)
        # distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceToOptStateBtwPops(popWithBrokenPrcs, pop1, 2, 2)
        # # Distance without broken PRC:
        # distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceBtwPops(popWithoutBrokenPrcs, pop1, 2, 1)
        # distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceBtwPops(popWithoutBrokenPrcs, pop1, 2, 2)
        # distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceToOptStateBtwPops(popWithoutBrokenPrcs, pop1, 2, 1)
        # distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1[j] = meanDistanceToOptStateBtwPops(popWithoutBrokenPrcs, pop1, 2, 2)
        # Broken PRC allows develstate to move closer to develstate in environment 1 when switch from environment 1 to environment 2, than when PRC is not broken

        ## Test when give totally new environment:
        # Find new environment to test that is different than both environment 1 and 2 evolved in by at least 50%:
        counts = 1
        differenceBtwEnvTwo = 0.0
        newEnv3State = zeros(Float64, 1, ENVS)
        while (differenceBtwEnvTwo <= 0.50) & (counts < MAXENVINPUT)
            EnvState2ToChange = randperm(ENVS)[1:convert(Int64,0.50*ENVS)]
            newEnv3State = copy(popWithoutBrokenPrcsAndTestNewEnv.founder.EnvState[:, 1])
            newEnv3State[EnvState2ToChange] = 1 .- popWithoutBrokenPrcsAndTestNewEnv.founder.EnvState[EnvState2ToChange, 1]
            differenceBtwEnvTwo = sum(abs.(newEnv3State - popWithoutBrokenPrcsAndTestNewEnv.founder.EnvState[:,2]))/ENVS
            counts += 1
        end
        # Save newEnv3State to print out for each independent trial:
        newEnv3StatesToPrint[:,j] = newEnv3State;
        # Test individuals in this new environment for stability:
        # Test for stability using previous develstate in given env 1 or 2 when introduce new env:
        # Env 1:
        map(x -> iterateindForStability(popWithoutBrokenPrcsAndTestNewEnv.individuals[x], 1, newEnv3State, 1, pop1.individuals[x]), 1:N)
        distanceWhenPolycombIntactInEnv1AndNewEnv[j] =  meanDistanceBtwPops(popWithoutBrokenPrcsAndTestNewEnv, pop1, 1, 1);
        percentPliantWhenPolycombIntactInEnv1AndNewEnv[j] =  phenotypicPliancyPercentage(popWithoutBrokenPrcsAndTestNewEnv, pop1, 1, 1);
        map(x -> iterateindForStability(popWithBrokenPrcsAndTestNewEnv.individuals[x], 1, newEnv3State, 1, pop1.individuals[x]), 1:N)
        distanceWhenPolycombBrokenInEnv1AndNewEnv[j] =  meanDistanceBtwPops(popWithBrokenPrcsAndTestNewEnv, pop1, 1, 1);
        percentPliantWhenPolycombBrokenInEnv1AndNewEnv[j] =  phenotypicPliancyPercentage(popWithBrokenPrcsAndTestNewEnv, pop1, 1, 1);
        # Save results:
        geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv1[:, (((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcsAndTestNewEnv,popWithoutBrokenPrcsAndTestNewEnv,1)
        # Number unstable individuals in Env 1:
        # When PRC intact:
        numStableIndWhenPrcIntactAndNewEnvReplaceEnv1[j] = length(findall(y -> y == true, map(x -> popWithoutBrokenPrcsAndTestNewEnv.individuals[x].stable[1], 1:N)))
        # When PRC broken:
        numStableIndWhenPrcBrokenAndNewEnvReplaceEnv1[j] = length(findall(y -> y == true, map(x -> popWithBrokenPrcsAndTestNewEnv.individuals[x].stable[1], 1:N)))
        # Env 2:
        map(x -> iterateindForStability(popWithoutBrokenPrcsAndTestNewEnv.individuals[x], 2, newEnv3State, 2, pop1.individuals[x]), 1:N)
        distanceWhenPolycombIntactInEnv2AndNewEnv[j] =  meanDistanceBtwPops(popWithoutBrokenPrcsAndTestNewEnv, pop1, 2, 2);
        percentPliantWhenPolycombIntactInEnv2AndNewEnv[j] =  phenotypicPliancyPercentage(popWithoutBrokenPrcsAndTestNewEnv, pop1, 2, 2);
        map(x -> iterateindForStability(popWithBrokenPrcsAndTestNewEnv.individuals[x], 2, newEnv3State, 2, pop1.individuals[x]), 1:N)
        distanceWhenPolycombBrokenInEnv2AndNewEnv[j] =  meanDistanceBtwPops(popWithBrokenPrcsAndTestNewEnv, pop1, 2, 2);
        percentPliantWhenPolycombBrokenInEnv2AndNewEnv[j] =  phenotypicPliancyPercentage(popWithBrokenPrcsAndTestNewEnv, pop1, 2, 2);
        # Save results:
        geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv2[:, (((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcsAndTestNewEnv,popWithoutBrokenPrcsAndTestNewEnv,2)
        # Number unstable individuals in Env 2:
        # When PRC intact:
        numStableIndWhenPrcIntactAndNewEnvReplaceEnv2[j] = length(findall(y -> y == true, map(x -> popWithoutBrokenPrcsAndTestNewEnv.individuals[x].stable[2], 1:N)))
        # When PRC broken:
        numStableIndWhenPrcBrokenAndNewEnvReplaceEnv2[j] = length(findall(y -> y == true, map(x -> popWithBrokenPrcsAndTestNewEnv.individuals[x].stable[2], 1:N)))
        ## Measurements when use initstate instead of develstate to start when switch envs:
        # Measure distance between individuals with PRC intact vs broken:
        # Env 1:
        # map(x -> iterateind2(popWithBrokenPrcsAndTestNewEnv.individuals[x], 1, newEnv3State), 1:N)
        # map(x -> iterateind2(popWithoutBrokenPrcsAndTestNewEnv.individuals[x], 1, newEnv3State), 1:N)
        # # Save results:
        # geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv1[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv1[:, (((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcsAndTestNewEnv,popWithoutBrokenPrcsAndTestNewEnv,1)
        # distanceWithNewEnvResultsReplaceEnv1[j] =  meanDistanceBtwPops(popWithBrokenPrcsAndTestNewEnv, popWithoutBrokenPrcsAndTestNewEnv, 1, 1)
        # ratioOfPcgBreakageToDistanceMoveInNewEnvResultsWhenReplaceEnv1[j] = meanRatioPcgBreakageToDistance(popWithBrokenPrcsAndTestNewEnv, popWithoutBrokenPrcsAndTestNewEnv, 1)
        # # Env 2:
        # map(x -> iterateind2(popWithBrokenPrcsAndTestNewEnv.individuals[x], 2, newEnv3State), 1:N)
        # map(x -> iterateind2(popWithoutBrokenPrcsAndTestNewEnv.individuals[x], 2, newEnv3State), 1:N)
        # Save results:
        # geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv2[:,(((j-1)*N)+1):(((j-1)*N)+N)], geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv2[:, (((j-1)*N)+1):(((j-1)*N)+N)] = gatherDevelStateResultsForBrokenAndIntact(popWithBrokenPrcsAndTestNewEnv,popWithoutBrokenPrcsAndTestNewEnv,2)
        # distanceWithNewEnvResultsReplaceEnv2[j] =  meanDistanceBtwPops(popWithBrokenPrcsAndTestNewEnv, popWithoutBrokenPrcsAndTestNewEnv, 2, 2)
        # ratioOfPcgBreakageToDistanceMoveInNewEnvResultsWhenReplaceEnv2[j] = meanRatioPcgBreakageToDistance(popWithBrokenPrcsAndTestNewEnv, popWithoutBrokenPrcsAndTestNewEnv, 2)

        # Measure distance between individuals when PRC intact and broken when in new env vs env evolved in:
        # For environment 1:
        # distanceWithNewEnvForEnv1WhenPrcIntactResults[j] =  meanDistanceBtwPops(pop1, popWithoutBrokenPrcsAndTestNewEnv, 1, 1)
        # distanceWithNewEnvForEnv1WhenPrcBrokenResults[j] =  meanDistanceBtwPops(pop1, popWithBrokenPrcsAndTestNewEnv, 1, 1)
        # # For env 2:
        # distanceWithNewEnvForEnv2WhenPrcIntactResults[j] =  meanDistanceBtwPops(pop1, popWithoutBrokenPrcsAndTestNewEnv, 2, 2)
        # distanceWithNewEnvForEnv2WhenPrcBrokenResults[j] =  meanDistanceBtwPops(pop1, popWithBrokenPrcsAndTestNewEnv, 2, 2)
        ## Test stability when only disrupt one PRC element:
        # popWith1BrokenPrcElement = deepcopy(pop1);
        # map(x -> breakpoly(popWith1BrokenPrcElement.individuals[x], 2), 1:N)
        # map(x -> iterateind2(popWith1BrokenPrcElement.individuals[x], 2, newEnv3State), 1:N)
        # numStableIndWhen1PrcElementBrokenAndNewEnvReplaceEnv2[j] = length(findall(y -> y == true, map(x -> popWith1BrokenPrcElement.individuals[x].stable[2], 1:N)))
        # # In environment 1:
        # map(x -> breakpoly(popWith1BrokenPrcElement.individuals[x], 1), 1:N)
        # map(x -> iterateind2(popWith1BrokenPrcElement.individuals[x], 1, newEnv3State), 1:N)
        # numStableIndWhen1PrcElementBrokenAndNewEnvReplaceEnv1[j] = length(findall(y -> y == true, map(x -> popWith1BrokenPrcElement.individuals[x].stable[1], 1:N)))
    end
    namesOfDataFrames = Vector{Symbol}(undef, N*INDTRIALS);
    namesOfDfsForSavingNewEnvs = Vector{Symbol}(undef, INDTRIALS);
    for i = 1:INDTRIALS
            namesOfDataFrames[((i-1)*N+1):(((i-1)*N)+N)] = Symbol.(Symbol(string("Trial",i,"-Ind")), 1:N)
            namesOfDfsForSavingNewEnvs[i] = Symbol(string("Trial",i))
    end
    savingGeneExpResults(geneExpAfterEvolEnv1, namesOfDataFrames, "geneExpAfterEvolEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpAfterEvolEnv2, namesOfDataFrames, "geneExpAfterEvolEnv2.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgIntactEnv1, namesOfDataFrames, "geneExpResultsPcgIntactEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgBrokenEnv1, namesOfDataFrames, "geneExpResultsPcgBrokenEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgIntactEnv2, namesOfDataFrames, "geneExpResultsPcgIntactEnv2.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgBrokenEnv2, namesOfDataFrames, "geneExpResultsPcgBrokenEnv2.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv1, namesOfDataFrames, "geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv1, namesOfDataFrames, "geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv2, namesOfDataFrames, "geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv2.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv2, namesOfDataFrames, "geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv2.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgIntactEnv2ReplaceEnv1, namesOfDataFrames, "geneExpResultsPcgIntactEnv2ReplaceEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgBrokenEnv2ReplaceEnv1, namesOfDataFrames, "geneExpResultsPcgBrokenEnv2ReplaceEnv1.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgIntactEnv1ReplaceEnv2, namesOfDataFrames, "geneExpResultsPcgIntactEnv1ReplaceEnv2.csv", folderForDataStorage)
    savingGeneExpResults(geneExpResultsPcgBrokenEnv1ReplaceEnv2, namesOfDataFrames, "geneExpResultsPcgBrokenEnv1ReplaceEnv2.csv", folderForDataStorage)
    # Save Polycomb mechanism control results:
    savingGeneExpResults(polycombStateAfterEvolEnv1, namesOfDataFrames, "polycombStateAfterEvolEnv1.csv", folderForDataStorage)
    savingGeneExpResults(polycombStateAfterEvolEnv2, namesOfDataFrames, "polycombStateAfterEvolEnv2.csv", folderForDataStorage)
    savingGeneExpResults(polycombVecAfterEvolEnv1, namesOfDataFrames, "polycombVecAfterEvolEnv1.csv", folderForDataStorage)
    savingGeneExpResults(polycombVecAfterEvolEnv2, namesOfDataFrames, "polycombVecAfterEvolEnv2.csv", folderForDataStorage)
    # Save new environment introduced (env 3) for each independent trial:
    savingGeneExpResults(newEnv3StatesToPrint, namesOfDfsForSavingNewEnvs, "newEnv3StatesForEachIndTrial.csv", folderForDataStorage)
    ## Results when PRC intact vs broken in evolved environments:
    outfile = joinpath(folderForDataStorage,"resultsWriteUp.txt")
    # writing to files is very similar:
    f = open(outfile, "w")
    println(f,"\nResults for independent trials = ", INDTRIALS, " and generations = ", GENS, " with randSeedNum = ", RANDSEEDNUM)
    println(f,"Evolved with ", PRCS, " PRC(S) and removed ", prcsToRemove)
    println(f,"\nFor when test for stability with develstate updated vs initstate:")
    println(f,"\nDistance when PcG broken & remain in Env 1 = ", meanOfResults(distanceWhenPolycombBrokenInEnv1))
    println(f,"\nDistance when PcG intact & remain in Env 1 = ", meanOfResults(distanceWhenPolycombIntactInEnv1))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenPolycombIntactInEnv1),removeNaNs(distanceWhenPolycombBrokenInEnv1))))
    println(f,"\nDistance when PcG broken and remain in Env 2 = ", meanOfResults(distanceWhenPolycombBrokenInEnv2))
    println(f,"\nDistance when PcG intact and remain in Env 2 = ", meanOfResults(distanceWhenPolycombIntactInEnv2))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenPolycombIntactInEnv2),removeNaNs(distanceWhenPolycombBrokenInEnv2))))
    println(f,"\nDistances when evolve in one env. but introduce other env.:")
    println(f,"\nDistance when evolve in env 1 and put in env 2 with PcG broken when compare to env 1 = ", meanOfResults(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1))
    println(f,"\nDistance when evolve in env 1 and put in env 2 with PcG intact when compare to env 1 = ", meanOfResults(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1),removeNaNs(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1))))
    println(f, "\nPercent phenotypically pliant for when evolve in env1 & put in env2 and compare to env1 with PcG broken = ", meanOfResults(percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1))
    println(f, "\nPercent phenotypically pliant for when evolve in env1 & put in env2 and compare to env1 with PcG intact = ", meanOfResults(percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1))
    println(f,"\npvalue % pliancy = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv1),removeNaNs(percentPliantWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv1))))
    println(f,"\nDistance when evolve in env 1 and put in env 2 with PcG broken when compare to env 2 = ", meanOfResults(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv2))
    println(f,"\nDistance when evolve in env 1 and put in env 2 with PcG intact when compare to env 2 = ", meanOfResults(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv2))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombIntactComparedToEnv2),removeNaNs(distanceWhenEvolveEnv1AndPutInEnv2WithPolycombBrokenComparedToEnv2))))
    println(f,"\nDistance when evolve in env 2 and put in env 1 with PcG broken when compare to env 2 = ", meanOfResults(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2))
    println(f,"\nDistance when evolve in env 2 and put in env 1 with PcG intact when compare to env 2 = ", meanOfResults(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2),removeNaNs(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2))))
    println(f, "\nPercent phenotypically pliant for when evolve in env2 & put in env1 and compare to env2 with PcG broken = ", meanOfResults(percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2))
    println(f, "\nPercent phenotypically pliant for when evolve in env2 & put in env1 and compare to env2 with PcG intact = ", meanOfResults(percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2))
    println(f,"\npvalue % pliancy = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv2),removeNaNs(percentPliantWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv2))))
    println(f,"\nDistance when evolve in env 2 and put in env 1 with PcG broken when compare to env 1 = ", meanOfResults(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv1))
    println(f,"\nDistance when evolve in env 2 and put in env 1 with PcG intact when compare to env 1 = ", meanOfResults(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv1))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombIntactComparedToEnv1),removeNaNs(distanceWhenEvolveEnv2AndPutInEnv1WithPolycombBrokenComparedToEnv1))))
    println(f,"\nResults when introduce new environment: ")
    println(f,"\nDistance when evovled in env 1 then PcG broken & put in new env vs env 1 develstate = ", meanOfResults(distanceWhenPolycombBrokenInEnv1AndNewEnv))
    println(f,"\nDistance when evovled in env 1 then PcG intact & put in new env vs env 1 develstate = ", meanOfResults(distanceWhenPolycombIntactInEnv1AndNewEnv))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenPolycombIntactInEnv1AndNewEnv),removeNaNs(distanceWhenPolycombBrokenInEnv1AndNewEnv))))
    println(f, "\nPercent phenotypically pliant for when evolve in env1 & put in new env vs env 1 develstate with PcG broken = ", meanOfResults(percentPliantWhenPolycombBrokenInEnv1AndNewEnv))
    println(f, "\nPercent phenotypically pliant for when evolve in env1 & put in new env vs env 1 develstate with PcG intact = ", meanOfResults(percentPliantWhenPolycombIntactInEnv1AndNewEnv))
    println(f,"\npvalue % pliancy = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(percentPliantWhenPolycombBrokenInEnv1AndNewEnv),removeNaNs(percentPliantWhenPolycombIntactInEnv1AndNewEnv))))
    println(f,"\nDistance when evovled in env 2 then PcG broken & put in new env vs env 2 develstate = ", meanOfResults(distanceWhenPolycombBrokenInEnv2AndNewEnv))
    println(f,"\nDistance when evovled in env 2 then PcG intact & put in new env vs env 2 develstate = ", meanOfResults(distanceWhenPolycombIntactInEnv2AndNewEnv))
    println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenPolycombIntactInEnv2AndNewEnv),removeNaNs(distanceWhenPolycombBrokenInEnv2AndNewEnv))))
    println(f,"\npvalue when broken PcG and left in env 1 vs when put in new env = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenPolycombBrokenInEnv1),removeNaNs(distanceWhenPolycombBrokenInEnv1AndNewEnv))))
    println(f,"\npvalue when broken PcG and left in env 2 vs when put in new env = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWhenPolycombBrokenInEnv2),removeNaNs(distanceWhenPolycombBrokenInEnv2AndNewEnv))))
    println(f, "\nPercent phenotypically pliant for when evolve in env2 & put in new env vs env 2 develstate with PcG broken = ", meanOfResults(percentPliantWhenPolycombBrokenInEnv2AndNewEnv))
    println(f, "\nPercent phenotypically pliant for when evolve in env2 & put in new env vs env 2 develstate with PcG intact = ", meanOfResults(percentPliantWhenPolycombIntactInEnv2AndNewEnv))
    println(f,"\npvalue % pliancy = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(percentPliantWhenPolycombBrokenInEnv2AndNewEnv),removeNaNs(percentPliantWhenPolycombIntactInEnv2AndNewEnv))))
    # When start at initstate instead of develstate:
    # println(f,"\nResults when PRC intact vs broken in evolved environments and start at initstate instead of develstate:")
    # println(f,"\nDistance to original develstate in env 1 when stay in env 1 & PRC broken = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenStayInEnv1))
    # println(f,"\nDistance to original develstate in env 1 when stay in env 1 & PRC intact = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv1WithIntactPrcWhenStayInEnv1))
    # println(f,"\nDistance to original develstate in env 2 when stay in env 2 & PRC broken = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenStayInEnv2))
    # println(f,"\nDistance to original develstate in env 2 when stay in env 2 & PRC intact = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv2WithIntactPrcWhenStayInEnv2))
    # println(f,"\nDistance to original develstate in env 1 when replace env 1 with env 2 & PRC broken = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\nDistance to original develstate in env 1 when replace env 1 with env 2 & PRC intact = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2),removeNaNs(distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2))))
    # println(f,"\nDistance to original develstate in env 2 when replace env 1 with env 2 & PRC broken = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\nDistance to original develstate in env 2 when replace env 1 with env 2 & PRC intact = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2),removeNaNs(distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2))))
    # println(f,"\nDistance to original develstate in env 1 when replace env 2 with env 1 & PRC broken = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\nDistance to original develstate in env 1 when replace env 2 with env 1 & PRC intact = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalDevelStateInGivenEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1),removeNaNs(distanceToOriginalDevelStateInGivenEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1))))
    # println(f,"\nDistance to original develstate in env 2 when replace env 2 with env 1 & PRC broken = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\nDistance to original develstate in env 2 when replace env 2 with env 1 & PRC intact = ", meanOfResults(distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalDevelStateInGivenEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1),removeNaNs(distanceToOriginalDevelStateInGivenEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1))))
    ## Results when PRC intact vs broken in evolved environments and distance from optstates:
    # In environment 1:
    # println(f,"\nDistance to optstate in env 1 when replace env 2 with env 1 & PRC broken = ", meanOfResults(distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\nDistance to optstate in env 1 when replace env 2 with env 1 & PRC intact = ", meanOfResults(distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv2WithEnv1),removeNaNs(distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv2WithEnv1))))
    # println(f,"\nDistance to optstate in env 2 when replace env 2 with env 1 & PRC broken = ", meanOfResults(distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\nDistance to optstate in env 2 when replace env 2 with env 1 & PRC intact = ", meanOfResults(distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv2WithEnv1),removeNaNs(distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv2WithEnv1))))
    # # In environment 2:
    # println(f,"\nDistance to optstate in env 1 when replace env 1 with env 2 & PRC broken = ", meanOfResults(distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\nDistance to optstate in env 1 when replace env 1 with env 2 & PRC intact = ", meanOfResults(distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalOptStateInEnv1WithBrokenPrcWhenReplaceEnv1WithEnv2),removeNaNs(distanceToOriginalOptStateInEnv1WithoutBrokenPrcWhenReplaceEnv1WithEnv2))))
    # println(f,"\nDistance to optstate in env 2 when replace env 1 with env 2 & PRC broken = ", meanOfResults(distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\nDistance to optstate in env 2 when replace env 1 with env 2 & PRC intact = ", meanOfResults(distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceToOriginalOptStateInEnv2WithBrokenPrcWhenReplaceEnv1WithEnv2),removeNaNs(distanceToOriginalOptStateInEnv2WithoutBrokenPrcWhenReplaceEnv1WithEnv2))))
    ## Results when introduce new environment:
    # distance btw intact and broken PRC when given new environment:
    # println(f,"\nResults when introduce new environment:")
    # println(f,"\nDistance btw individuals with PRC intact vs broken when given new env. in place of env 1 = ", meanOfResults(distanceWithNewEnvResultsReplaceEnv1))
    # println(f,"\nDistance btw individuals with PRC intact vs broken when given new env. in place of env 2 = ", meanOfResults(distanceWithNewEnvResultsReplaceEnv2))
    # println(f,"\nRatio of PcG genes suppressed/ENVS to distance btw intact & broken when given new env. in place of env 1 = ",meanOfResults(ratioOfPcgBreakageToDistanceMoveInNewEnvResultsWhenReplaceEnv1))
    # println(f,"\nRatio of PcG genes suppressed/ENVS to distance btw intact & broken when given new env. in place of env 2 = ",meanOfResults(ratioOfPcgBreakageToDistanceMoveInNewEnvResultsWhenReplaceEnv2))
    # println(f,"\nResults when introduce new environment when PRC intact vs broken:")
    # println(f,"\nDistance btw individuals with PRC broken in new env vs env 1 = ", meanOfResults(distanceWithNewEnvForEnv1WhenPrcBrokenResults))
    # println(f,"\nDistance btw individuals with PRC intact in new env vs env 1 = ", meanOfResults(distanceWithNewEnvForEnv1WhenPrcIntactResults))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWithNewEnvForEnv1WhenPrcIntactResults),removeNaNs(distanceWithNewEnvForEnv1WhenPrcBrokenResults))))
    # println(f,"\nDistance btw individuals with PRC broken in new env vs env 2 = ", meanOfResults(distanceWithNewEnvForEnv2WhenPrcBrokenResults))
    # println(f,"\nDistance btw individuals with PRC intact in new env vs env 2 = ", meanOfResults(distanceWithNewEnvForEnv2WhenPrcIntactResults))
    # println(f,"\npvalue = ", pvalue(ApproximateTwoSampleKSTest(removeNaNs(distanceWithNewEnvForEnv2WhenPrcIntactResults),removeNaNs(distanceWithNewEnvForEnv2WhenPrcBrokenResults))))
    ## Stability results:
    println(f,"\nStability results:")
    println(f,"\nNumber stable individuals with PRC broken when given new env. in place of env 1 = ", mean(numStableIndWhenPrcBrokenAndNewEnvReplaceEnv1))
    println(f,"\nNumber stable individuals with PRC intact when given new env. in place of env 1 = ", mean(numStableIndWhenPrcIntactAndNewEnvReplaceEnv1))
    println(f,"\nNumber stable individuals with PRC broken when given new env. in place of env 2 = ", mean(numStableIndWhenPrcBrokenAndNewEnvReplaceEnv2))
    println(f,"\nNumber stable individuals with PRC intact when given new env. in place of env 2 = ", mean(numStableIndWhenPrcIntactAndNewEnvReplaceEnv2))
    # print("\nNumber stable individuals with only 1 PRC target gene broken when given new env. in place of env 1 = ", mean(numStableIndWhen1PrcElementBrokenAndNewEnvReplaceEnv1))
    # print("\nNumber stable individuals with only 1 PRC target gene broken when given new env. in place of env 2 = ", mean(numStableIndWhen1PrcElementBrokenAndNewEnvReplaceEnv2))
    ## Fitness results:
    println(f,"\nFitness of pop. in [env1; env2] before evolution = ", mean(fitnessOfIndividualsAtStartOfEvolution, dims = 2))
    println(f,"\nFitness of pop. in [env1; env2] after final generation = ", mean(fitnessOfIndividualsAtEndOfEvolution, dims = 2))
    println(f,"\nOverall fitness of pop. before & after evolution = ", mean(overallFitnessBeforeEvolution), " & ", mean(overallFitnessAfterEvolution))
    ## Genes suppressed:
    println(f,"\nAvg genes suppressed by PRCs = ", mean(averageNumOfGenesSuppressedByPrcs, dims=2))
    close(f)
end







#### Test what happens to resulting develstate when introduce different initstates:
#Pkg.build("GR") # if restarting julia
using Plots
randSeedNumber = 171#rand(1:1000);
println("random seed = ", randSeedNumber)
Random.seed!(randSeedNumber);
pop1 = genpop();
println("fitnessOfIndividualsAtStartOfEvolution = ", mean(map(x -> pop1.individuals[x].fitnessUnderEachEnv, 1:N)))
println("overallFitnessBeforeEvolution = ", mean(map(x -> pop1.individuals[x].fitness, 1:N)))
for i = 1:GENS
    update(pop1);
end
println("fitnessOfIndividualsAtEndOfEvolution = ", mean(map(x -> pop1.individuals[x].fitnessUnderEachEnv, 1:N)))
println("overallFitnessAfterEvolution = ", mean(map(x -> pop1.individuals[x].fitness, 1:N)))
# Find polycomb mechanism results after evolution:
println("averageNumOfGenesSuppressedByPrcs = ", mean(map(z -> map(y -> length(findall(x -> x == 0, pop1.individuals[z].polycombstate[:,y])), 1:2), 1:N)))
distanceForFitness = map(x->sum(abs.(pop1.individuals[x].develstate[:,2] - pop1.individuals[x].optstate[:,2]))/G,1:N)
println("\nMean distance btw develstate & optstate of env 2 for original initstate = ", round(mean(distanceForFitness),digits=3), "\nstd of distance btw develstate & optstate of env 2 for original initstate = ", round(std(distanceForFitness),digits=3))
### For initstate testing:
meanDistanceForFitnessNewInitState = zeros(1000)
stdDistanceForFitnessNewInitState = zeros(1000)
for i = 1:1000
    newInitState = rand(0:1,G)
    map(x->pop1.individuals[x].initstate = newInitState, 1:N)
    map(x->iterateind2(pop1.individuals[x], 2, pop1.founder.EnvState[:,2]), 1:N)
    distanceForFitnessNewInitState = map(x->sum(abs.(pop1.individuals[x].develstate[:,2] - pop1.individuals[x].optstate[:,2]))/G,1:N)
    meanDistanceForFitnessNewInitState[i] = mean(distanceForFitnessNewInitState)
    stdDistanceForFitnessNewInitState[i] = std(distanceForFitnessNewInitState)
    #println("\n Mean distance btw develstate & optstate of env 2 for new initstate = ", mean(distanceForFitnessNewInitState), "\nstd of distance btw develstate & optstate of env 2 for new initstate = ", std(distanceForFitnessNewInitState))
end
plot(histogram(meanDistanceForFitnessNewInitState, title="mean"),histogram(stdDistanceForFitnessNewInitState, title="std"),layout=(2,1))
println("\nMean distance btw develstate & optstate of env 2 for original initstate = ", round(mean(distanceForFitness),digits=3), "\nstd of distance btw develstate & optstate of env 2 for original initstate = ", round(std(distanceForFitness),digits=3))
println("Range of means for when vary initstate = ", round(minimum(meanDistanceForFitnessNewInitState),digits=3), " to ", round(maximum(meanDistanceForFitnessNewInitState),digits=3))
println("Range of stds for when vary initstate = ", round(minimum(stdDistanceForFitnessNewInitState),digits=3), " to ", round(maximum(stdDistanceForFitnessNewInitState),digits=3))
### For environment testing of variability of develstate output:
map(x->pop1.individuals[x].initstate = copy(pop1.founder.initstate), 1:N)
meanDistanceForFitnessNewEnv = zeros(1000)
stdDistanceForFitnessNewEnv = zeros(1000)
for i = 1:1000
    counts = 1
    differenceBtwEnvTwo = 0.0
    newEnv3State = zeros(Float64, 1, ENVS)
    while (differenceBtwEnvTwo <= 0.50) & (counts < MAXENVINPUT)
        EnvState2ToChange = randperm(ENVS)[1:convert(Int64,0.50*ENVS)]
        newEnv3State = copy(pop1.founder.EnvState[:, 1])
        newEnv3State[EnvState2ToChange] = 1 .- pop1.founder.EnvState[EnvState2ToChange, 1]
        differenceBtwEnvTwo = sum(abs.(newEnv3State - pop1.founder.EnvState[:,2]))/ENVS
        counts += 1
    end
    map(x->iterateind2(pop1.individuals[x], 2, newEnv3State), 1:N)
    distanceForFitnessNewEnv = map(x->sum(abs.(pop1.individuals[x].develstate[:,2] - pop1.individuals[x].optstate[:,2]))/G,1:N)
    meanDistanceForFitnessNewEnv[i] = mean(distanceForFitnessNewEnv)
    stdDistanceForFitnessNewEnv[i] = std(distanceForFitnessNewEnv)
end
plot(histogram(meanDistanceForFitnessNewEnv, title="mean"),histogram(stdDistanceForFitnessNewEnv, title="std"),layout=(2,1))
println("Range of means for when vary env = ", round(minimum(meanDistanceForFitnessNewEnv),digits=3), " to ", round(maximum(meanDistanceForFitnessNewEnv),digits=3))
println("Range of stds for when vary env = ", round(minimum(stdDistanceForFitnessNewEnv),digits=3), " to ", round(maximum(stdDistanceForFitnessNewEnv),digits=3))

println("\n p-value for means of new inistates vs env states = ",pvalue(ApproximateTwoSampleKSTest(meanDistanceForFitnessNewEnv,meanDistanceForFitnessNewInitState)))
println("\n p-value for stds of new inistates vs env states = ",pvalue(ApproximateTwoSampleKSTest(stdDistanceForFitnessNewEnv,stdDistanceForFitnessNewInitState)))
