function genmeasure()
    time = zeros(Int64,NUMMEAS) # create initial vectors with length equal to the total number of times going to make measurements of population properties
    fitness = zeros(Float64,NUMMEAS)
    fitnessstd = zeros(Float64,NUMMEAS)
    fitnessUnderEachEnv = zeros(Float64,NUMMEAS,DIFENVS)
    fitnessUnderEachEnvStd = zeros(Float64,NUMMEAS,DIFENVS)
    robustness = zeros(Float64,NUMMEAS,DIFENVS)
    robustnessstd = zeros(Float64,NUMMEAS,DIFENVS)
    envRobustness = zeros(Float64,NUMMEAS,DIFENVS)
    envRobustnessStd = zeros(Float64,NUMMEAS,DIFENVS)
    pathlength = zeros(Float64,NUMMEAS,DIFENVS)
    pathlengthstd = zeros(Float64,NUMMEAS,DIFENVS)
    indtypes = zeros(Float64,NUMMEAS)
    inittypes = zeros(Float64,NUMMEAS)
    develtypes = zeros(Float64,NUMMEAS,DIFENVS)
    pcgVecTypes = zeros(Float64, NUMMEAS)
    pcgStateTypes = zeros(Float64,NUMMEAS,DIFENVS)
    opttypes = zeros(Float64,NUMMEAS,DIFENVS)
    connectivity = zeros(Float64, NUMMEAS)

    # set Measure type with the initial values above
    Measure(time, fitness, fitnessstd, fitnessUnderEachEnv, fitnessUnderEachEnvStd, robustness, robustnessstd, envRobustness,
            envRobustnessStd, pathlength, pathlengthstd, indtypes, inittypes,
            develtypes, pcgVecTypes, pcgStateTypes, opttypes, connectivity)
end

# Set measurement functions depending on type of measurement wanting to take:
measureFuncs = Dict{String, Function}();
measureFuncs["time"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> [t];
measureFuncs["indtypes"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> [length(unique(map(x->x.network, pop.individuals)))];
measureFuncs["inittypes"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> [length(unique(map(x->x.initstate, pop.individuals)))];
measureFuncs["connectivity"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> [mean(map(x->count(z->(z==0),x),map(x->x.network, pop.individuals)))];
measureFuncs["develtypes"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> map(y -> length(unique(map(x->x.develstate[:,y], pop.individuals))), 1:DIFENVS) # old code was: mean(map(x->sum(abs.(pop.founder.optstate - x.develstate),dims=1)/G, pop.individuals)) # mean of distance from optimum to each individual's development vector (S hat)
measureFuncs["pcgVecTypes"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> begin pcgSusceptible = map(x->count(y->y==1,sum(x.polycombvec, dims = 2)), pop.individuals); # finding how many genes per individual is susceptible to at least one PRC
                                                        if length(findall(x -> x == 1, pcgSusceptible)) == 0
                                                            [0]
                                                        else
                                                            [mean(pcgSusceptible[findall(x -> x == 1, pcgSusceptible)])] # mean of number of susceptible genes for all individuals -> find all individuals that have at least 1 susceptible gene
                                                        end
                                                    end
measureFuncs["pcgStateTypes"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> begin pcgSupressed = map(x->sum(z->z==0, x.polycombstate, dims=1), pop.individuals);
                                                        pcgActuallySupressed = map(y->map(x->pcgSupressed[x][y],1:size(pcgSupressed)[1]),1:DIFENVS);
                                                        numGenesSuppressed = zeros(Float64,DIFENVS);
                                                        for i = 1:DIFENVS
                                                            if length(findall(x -> x > 0, pcgActuallySupressed[i])) == 0 # find if at least 1 gene is suppressed by polycomb
                                                                numGenesSuppressed[i] = 0.
                                                            else
                                                                numGenesSuppressed[i] = mean(pcgActuallySupressed[i][findall(x -> x > 0, pcgActuallySupressed[i])])
                                                            end
                                                        end
                                                        numGenesSuppressed
                                                    end # mean of number of genes supressed by polycomb (not taking into account genes that are not suppressed by polycomb)
measureFuncs["opttypes"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> map(y -> length(unique(map(x->x.optstate[:,y], pop.individuals))), 1:DIFENVS)
measureFuncs["fitness"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> [mean(map(x->x.fitness, pop.individuals))]
measureFuncs["fitnessstd"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> [std(map(x->x.fitness, pop.individuals))]
measureFuncs["pathlength"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> mean(map(x->x.pathlength, pop.individuals))
measureFuncs["pathlengthstd"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> std(map(x->x.pathlength, pop.individuals))
measureFuncs["robustness"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> mean(robustvect)
measureFuncs["robustnessstd"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> std(robustvect)
measureFuncs["envRobustness"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> mean(envrobustvect)
measureFuncs["envRobustnessStd"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> std(envrobustvect)
measureFuncs["fitnessUnderEachEnv"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> mean(map(x->x.fitnessUnderEachEnv, pop.individuals))
measureFuncs["fitnessUnderEachEnvStd"] = (pop::Population, t::Int64, robustvect::Array{Array{Float64,1},1}, envrobustvect::Array{Array{Float64,1},1}) -> std(map(x->x.fitnessUnderEachEnv, pop.individuals))

function measure(pop::Population, m::Measure, t::Int64, n::Int64, mutatedEnvMatrix::Array{Array{Int64,2},1}, founderGeneticRobustness::Vector{Float64},
    founderEnvironmentalRobustness::Vector{Float64}, vecMeasurementsToTake::Array{String,1}, measureTypes::Dict{String,Array})
# Take in a Population type, pop; Measure Type, m;
# Time, t, representing the current generation (so step width of 1, like 1:1:1000);
# and the number of the particular measurement being taken, n (1:total # of
# measurements to take (TOTALMEAS)):
    if "envRobustness" in vecMeasurementsToTake
        map(i -> environmentalRobustness(i, mutatedEnvMatrix, founderEnvironmentalRobustness), pop.individuals)
        envrobustvect = map(i -> pop.individuals[i].envRobustness, 1:N)
    else
        envrobustvect = zeros(Float64,N,DIFENVS)
        envrobustvect = map(i -> envrobustvect[i,:], 1:N)
    end
    if "robustness" in vecMeasurementsToTake
        map(i -> geneticRobustness(i, founderGeneticRobustness), pop.individuals)
        robustvect = map(i -> pop.individuals[i].robustness, 1:N)
    else
        robustvect = zeros(Float64,N,DIFENVS)
        robustvect = map(i -> robustvect[i,:], 1:N)
    end
    # calculate measurements that user inputs to calculate:
    for i=1:length(vecMeasurementsToTake)
        measureTypes[vecMeasurementsToTake[i]][n,:] = measureFuncs[vecMeasurementsToTake[i]](pop, t, robustvect, envrobustvect);
    end
    # store measurements in Measure:
    m.time[n] = measureTypes["time"][n,:][1]
    m.fitness[n] = measureTypes["fitness"][n,:][1]
    m.fitnessstd[n] = measureTypes["fitnessstd"][n,:][1]
    m.fitnessUnderEachEnv[n,:] = measureTypes["fitnessUnderEachEnv"][n,:]
    m.fitnessUnderEachEnvStd[n,:] = measureTypes["fitnessUnderEachEnvStd"][n,:]
    m.robustness[n,:] = measureTypes["robustness"][n,:]
    m.robustnessstd[n,:] = measureTypes["robustnessstd"][n,:]
    m.envRobustness[n,:] = measureTypes["envRobustness"][n,:]
    m.envRobustnessStd[n,:] = measureTypes["envRobustnessStd"][n,:]
    m.pathlength[n,:] = measureTypes["pathlength"][n,:]
    m.pathlengthstd[n,:] = measureTypes["pathlengthstd"][n,:]
    m.indtypes[n] = measureTypes["indtypes"][n,:][1]
    m.inittypes[n] = measureTypes["inittypes"][n,:][1]
    m.develtypes[n,:] = measureTypes["develtypes"][n,:]
    m.pcgVecTypes[n] = measureTypes["pcgVecTypes"][n,:][1]
    m.pcgStateTypes[n,:] = measureTypes["pcgStateTypes"][n,:]
    m.opttypes[n,:] = measureTypes["opttypes"][n,:]
    m.connectivity[n] = measureTypes["connectivity"][n,:][1]
end


function measureUsingSaurabhCcodeVersionForEnvRobust(pop::Population, m::Measure, t::Time, n::Int64, mutatedEnvMatrix::Array{Float64,2}, founderGeneticRobustness::Float64,
    founderEnvironmentalRobustness::Float64, vecMeasurementsToTake::Array{String,1}, measureTypes::Dict{String,Array{T,1} where T})
# Take in a Population type, pop; Measure Type, m;
# Time, t, representing the current generation (so step width of 1, like 1:1:1000);
# and the number of the particular measurement being taken, n (1:total # of
# measurements to take (TOTALMEAS)):
    if "envRobustness" in vecMeasurementsToTake
        envrobustvect = pmap(i -> environmentalRobustnessSaurabhCcodeVersion(i, mutatedEnvMatrix, founderEnvironmentalRobustness, pop.founder), pop.individuals)
        for i = 1:N # loop to set robustness of individuals **Is this NOT done else where in individuals.jl file?? -ML 3/26/19
            pop.individuals[i].envRobustness = envrobustvect[i]
        end
    else
        envrobustvect = zeros(Float64,N)
    end
    if "robustness" in vecMeasurementsToTake
        robustvect = pmap(i -> geneticRobustness(i, founderGeneticRobustness), pop.individuals)
        for i = 1:N # loop to set robustness of individuals **Is this NOT done else where in individuals.jl file?? -ML 3/26/19
            pop.individuals[i].robustness = robustvect[i]
        end
    else
        robustvect = zeros(Float64,N)
    end
    # calculate measurements that user inputs to calculate:
    for i=1:length(vecMeasurementsToTake)
        measureTypes[vecMeasurementsToTake[i]][n] = measureFuncs[vecMeasurementsToTake[i]](pop, t, robustvect, envrobustvect);
    end
    # store measurements in Measure:
    m.time[n] = measureTypes["time"][n]
    m.fitness[n] = measureTypes["fitness"][n]
    m.fitnessstd[n] = measureTypes["fitnessstd"][n]
    m.robustness[n] = measureTypes["robustness"][n]
    m.robustnessstd[n] = measureTypes["robustnessstd"][n]
    m.envRobustness[n] = measureTypes["envRobustness"][n]
    m.envRobustnessStd[n] = measureTypes["envRobustnessStd"][n]
    m.pathlength[n] = measureTypes["pathlength"][n]
    m.pathlengthstd[n] = measureTypes["pathlengthstd"][n]
    m.indtypes[n] = measureTypes["indtypes"][n]
    m.inittypes[n] = measureTypes["inittypes"][n]
    m.develtypes[n] = measureTypes["develtypes"][n]
    m.pcgVecTypes[n] = measureTypes["pcgVecTypes"][n]
    m.pcgStateTypes[n] = measureTypes["pcgStateTypes"][n]
    m.opttypes[n] = measureTypes["opttypes"][n]
    m.connectivity[n] = measureTypes["connectivity"][n]
end


function standardErrorOfMean(x::Vector)
# Input x is vector of some sort of values,
# like mean environmental robustness vector
# for all independent trials. Then this function
# returns standard error of the mean of the
# given vector, x.
    return(std(x)/sqrt(length(x)))
end
