# ---------------------------
# CONSTRUCTORS
# ---------------------------
function Individual(network::Matrix{Float64}, initstate::Vector{Float64})
    Individual(copy(network), copy(initstate), zeros(Float64,(G,ENVS)), zeros(Float64,ENVS,DIFENVS),
                zeros(Float64,G,DIFENVS), zeros(Float64,G,DIFENVS), falses(DIFENVS), 0., zeros(Float64,DIFENVS), zeros(Int64,DIFENVS),
                    zeros(Float64,G,PRCS), ones(Float64,G,DIFENVS), zeros(Float64,DIFENVS), zeros(Float64,DIFENVS))
    # Individual composed of:
    # Individual(network, initstate, EnvIndInter, EnvState,
    #         develstate, optstate, stable, fitness, fitnessUnderEachEnv, pathlength,
    #                 polycombvec, polycombstate, robustness, envRobustness)
end

function Individual(network::Matrix{Float64})
    Individual(copy(network), zeros(Float64,G))
end

function Individual()
    Individual(zeros(G,G))
end


# ----------------------------------------------
# METHODS: Generate Individuals
# ----------------------------------------------
function geninds(founder::Individual)
# If given an individual to use as the founder:
# Initialize a Population by creating
# N individuals that are all a copy of founder
    indvect1 = Array{Individual{Float64}}(undef, N) # Julia assigns Individual type as whole as Float64
    for i = 1:N
        if RANDPOP # generate random population that is NOT a copy of the founder
            indvect1[i] = stableind(founder.initstate)
        else # All individuals in population are copy of a founder
            indvect1[i] = deepcopy(founder)
        end
    end
    return indvect1
end


function geninds()
# If not given an individual to use as the founder:
# Generate N individuals by generating a
# stable founder first and then copy this
# founder N times
    founder = genfounder()

    indvect1 = Array{Individual{Float64}}(undef, N) # Julia assigns Individual type as whole as Float64
    for i = 1:N
        if RANDPOP # generate random population that is NOT a copy of the founder
            indvect1[i] = stableind(founder.initstate)
        else # All individuals in population are copy of a founder
            indvect1[i] = deepcopy(founder)
        end
    end
    return indvect1
end

# ----------------------------------------------
# METHODS: Generate Founder Step
# ----------------------------------------------
function genfounder()
# generate a founding individual whose
# developmentally stable state is equivalent to
# its optimal state (however rand
# assumes user not putting in optimal state vector)
# and thereby determines the optimal state for the population,
# i.e. the optimal state is set as the final state of founder
# once stability is reached within MAXCONV iterations
    if isempty(INP) & isempty(OPT)# **this is what currently runs as of 7/25/18 bc INP is NOT set in constants now (bc caused problems when trying to loop through to generate many different populations with different INP1 because had to include("../input/constants.jl") everytime) - ML 7/25/18
        founder = randstableind()
        # print(" initstate & optstate empty")
    elseif isempty(OPT) #this what used to run bc INP was set in constants (before 7/25/18)
        founder = stableind(INP) # founder should be stable at this point
        # print(" inistate NOT empty")
    elseif isempty(INP)
        # print(" optstate NOT empty")
        founder = randstableind()
    else
        # print(" initstate & optstate NOT empty")
        founder = stableind(INP)
    end
    return founder
end

function randstableind() # set Sw vector
# generate a random Individual that is
# developmentally stable by initializing the state vector as rand(0:1,G)*2-1 and passing through stableind function
# initialzing the initial state vector and generating a stable individual
    randind = stableind(convert(Vector{Float64},rand(0:1,G))) # assumes user not inputting optstate into stableind and then generates initial state
    return randind
end

function stableind(initstate::Vector{Float64}) # used if only do not input optstate into function; optstate will be empty vector
    if isempty(OPT) & isempty(WMAT) #this is what usually runs as of Aug 2019 because W & optstate are empty in constants
        ind = stableind(initstate,(Float64)[],Array{Float64,2}(undef,0,0))
        # print(" wmat & optstate NOT given")
    elseif isempty(OPT)
        ind = stableind(initstate,(Float64)[],WMAT)
        # print(" optstate NOT given, w mat given")
    elseif isempty(WMAT)
        ind = stableind(initstate,OPT,Array{Float64,2}(undef,0,0))
        # print(" wmat NOT given, optstate given")
    else
        ind = stableind(initstate,OPT,WMAT)
        # print(" given optstate & wmatrix")
    end
    return ind
end


function stableind(initstate::Vector{Float64},optstate::Vector{Float64},wmatrix::Matrix{Float64})
# This function is used to generate a founder that is stable:
# generate a random Individual that is developmentally stable by
# giving an initial individual state vector and can optionally give
# optimal state vector if desire. This function generates an individual's
# Gene Regulatory Network (W), Environment State vector (Se) and
# Environment Gene interaction matrix (W x E).
    if isempty(wmatrix)
        randind = Individual() # generating object of type individual
        randind.initstate = copy(initstate)
        # print(" wmat NOT given in stableind(inistate,optstate,wmat) function")
    else
        # print(" wmat given in stableind(inistate,optstate,wmat) function")
        randind = Individual() # generating object of type individual
        randind.initstate = copy(initstate)
        randind.network = copy(wmatrix)
    end
    # what happens if this wmatrix and initstate isn't stable for saurabh version but for extended?? 3/21/19 -ML
    # so maybe always want to generate Saurabh population first and use this Wmatrix to find stable extended population founder
    # b/c can alter WxE interaction matrix and/or Se vector when searching for stability with given initstate vector and W matrix
    countsForInitstates = 1
    while (randind.stable!=trues(DIFENVS)) & (countsForInitstates < MAXINITSTATES)
        while ((randind.stable[1]!=true) | (if isempty(optstate); false; else; randind.develstate!=optstate; end))
            # right in OR statment is only for when user inputs an optstate,
            # but in PcG model optstate is always empty because set as such in constants
            if isempty(wmatrix) # if not W matrix (network) given, then generate random one
                randind.network = zeros(Float64,G,G)
                if IPANETWORKS == false
                    for i=1:G^2
                        if rand()<C # connectivity probability for individuals' W (Gaussian matrices)
                            if INDWEIGHTS == "discrete" # individual weights; if "discrete" then sample {1, -1, 0} via rand(-1:1)
                                #randind.network[i] = [-1,1][rand(1:2)]
                                randind.network[i] = rand(-1:1)
                            elseif INDWEIGHTS == "gaussian" # if "gaussian" then sample -inf to inf (is this correct??? 2/28/18 -ML --> yes) via randn()
                                randind.network[i] = randn()
                            else
                                error(string("wrong specification of individual weights:",
                                       " use discrete or gaussian"))
                            end
                        end
                    end
                else
                    randind.network = generateWmatUsingIpaNetwork(IPAFILENAME)
                end
            end

            # Initialize Environment Gene Network Interaction matrix (W x E)
            randind.EnvIndInter = zeros(Float64, G, ENVS)
            if !SAURABVERSION # if running Saurab's older version of PcG, then do not have environmentalStateVec & environmentInteractionMatrix (thus equal to zero, which are set in Individual function above)
                for i=1:G
                    if rand()<CENMAT # proportion of environments that can possibly affect a gene(s)
                        for j=1:ENVS
                            if rand()<CENGENE # proportion of genes that can be affected by a given environment
                                randind.EnvIndInter[i, j] = randn() # effect of environment j on gene i
                            end
                        end
                    end
                end
            end

            # Initialize Environment State vector (environment present or not)
            if SAURABVERSION # if running Saurab's older version of PcG, then do not have environmentalStateVec (thus equal to zero)
                randind.EnvState = zeros(Float64,ENVS,DIFENVS)
            else
                envState1 = convert(Vector{Float64},rand(0:1,ENVS))
                envStatesVectors = map(x -> envState1, 1:DIFENVS)
                randind.EnvState = hcat(envStatesVectors...)
            end
            iterateindfound(randind, 1)
        end

        # For different environmental states:
        # First find other environmental inputs (EnvState) that are stable for randind
        # but 70% different than envState of randind.envState[:,1]
        # assuming that initstate are the same for both environmental state cases
        if DIFENVS > 1
            for i = 2:DIFENVS
                counts = 1
                while (randind.stable[i] == false) & (counts < MAXENVINPUT)
                    EnvState2ToChange = randperm(ENVS)[1:convert(Int64,PERCENTDIFENVS*ENVS)]
                    randind.EnvState[EnvState2ToChange, i] = 1 .- randind.EnvState[EnvState2ToChange, 1]
                    iterateindfound(randind, i) # develstate, stable, pathlength, and polycombstate get set in this function for this particular environmental state, i
                    counts += 1
                end
            end
        end
        countsForInitstates += 1
    end
    if countsForInitstates >= MAXINITSTATES
        randind = stableind(convert(Vector{Float64},rand(0:1,G)))
    end
    # set optstate
    if isempty(optstate)
        randind.optstate[:,1] = copy(randind.develstate[:,1]) # set optstate of first environment so this environment is one that is already optimum
    else # if user inputs optstate then develstate and optstate will not be equal
        randind.optstate = copy(optstate)
    end
    # Find and set optstate(s) in other environmental states that is at least 40% different than founder opstate in environment state 1
    if DIFENVS > 1
        for j = 2:DIFENVS
            distanceOptStates = 0.0
            countNum = 1
            while (distanceOptStates < PERCENTDIFOPTSTATES) & (countNum < MAXOPTSTATES)
                randind.optstate[:,j] = copy(randind.optstate[:,1])
                optState2ToChange = randperm(G)[1:convert(Int64,PERCENTDIFGENES*G)]
                randind.optstate[optState2ToChange,j] = abs.(rand() .- randind.optstate[optState2ToChange,1])
                distanceOptStates = sum(abs.(randind.optstate[:,1] - randind.optstate[:,j]))/G
                countNum += 1
            end
            if countNum > MAXOPTSTATES
                randind = stableind(initstate)
            end
        end
    end
    # calculate fitness
    fitnesseval(randind)
    fitnessUnderEachEnvEval(randind)
    return randind
end


function iterateindfound(me::Individual, envState::Int64)
# iterate founder individual from its initial state to
# their developmental state and store the results
# in the Individual object under only the first environmental state
# This function used to find a developmentally stable founder to use for initial population
    currstate = copy(me.initstate)
    vectAD = zeros(Float64,MAXCONV)
    for i=1:MAXCONV
        stateupdate = me.network*currstate
        stateupdateenv = me.EnvIndInter*me.EnvState[:,envState]
        stateupdate = stateupdate + stateupdateenv
        stateupdate = 1.0./(1.0 .+ exp.(-SIGSTR.*stateupdate))
        tempdiff = currstate - stateupdate
        actualDist = sum(abs.(tempdiff))/G

        #sliding window - added 7/5/18 - SG
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.0))
            stableBool = slidingWindow(i,vectAD)
            if distBool && stableBool
                me.stable[envState] = true
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = copy(i) #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable[envState] = false
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = copy(MAXCONV) #changed to copy() - SG 6/21/18
            end
        end
        currstate = copy(stateupdate)
    end
end


function slidingWindow(t::Int64, x::Array{Float64,1})
# Measures if individual has reached stability using a sliding window of t-(TAU-1) to t and
# taking the normalized Hamming distance of the measurements recorded and tests if they
# fulfill the criteria of stability - i.e. the average Hamming distances in a sliding
# window of length TAU is less than epsilon, which equals 10^-4.
# Written by Sarah Garaff July 5, 2018
    distances = 0.0
    for i = (t-(TAU-1)):t
        distances += x[i] #add up each of the distances between measurements
    end
    std = distances/TAU #find std by divding by TAU
    if std <= EPSI # stable if normalized variance is less than epsilon (EPSI), which equals 10^(-6)
         return true # it's stable
    else
        return false
    end
end


# -----------------------------------------
# METHODS: Generate Random Environments
# to measure Environmental Robustness Step
# -----------------------------------------
function generateMutatedEnvsMatrix(me::Individual, envState::Int64)
# If given an individual type (founder individual in main run julia file),
# then generate matrix with all possible environment
# state vectors (extended version) OR initial state vectors (Saurabh version)
# that are only different from initial environment state vector (extended
# version) or initial state vector (Saurabh version) of given individual
# by ALTERENV percentage of environments being changed.
# If running Saurab's version, then change the initstate vector of the
# given individual by ALTERENV*G (# genes).
# If running extended version, then change EnvState of the given
# individual by ALTERENV*ENVS (# of environments).
    if SAURABVERSION
        mutatedInitStates = collect(combinations(collect(1:length(me.initstate)), convert(Int64, ALTERENV*G)))
        mutatedInitStateIndx = sample(1:length(mutatedInitStates), ENVCHANGES, replace = false)
        mutatedInitStateMatrix = zeros(Float64, G, ENVCHANGES)
        for i=1:ENVCHANGES
            copyInitState = copy(me.initstate) # copy founder.initstate without changing founder.initstate or me.initstate
            copyInitState[mutatedInitStates[mutatedInitStateIndx[i]]] = 1 .- copyInitState[mutatedInitStates[mutatedInitStateIndx[i]]]
            mutatedInitStateMatrix[:,i] = copyInitState
        end
        return(mutatedInitStateMatrix)
    else
        mutatedEnvs = collect(combinations(collect(1:length(me.EnvState[:,envState])), convert(Int64, ALTERENV*ENVS))) # all combinations of
            # changing only 20% (ALTERENV) of the elements of Environment state vector (ENVS)
        mutatedEnvIndx = sample(1:length(mutatedEnvs), ENVCHANGES, replace = false) # find indexes of environments want to mutate (mutatedEnvs)
            # for the number of times want to change environment (ENVCHANGES)
        mutatedEnvMatrix = zeros(Float64, ENVS, ENVCHANGES)
        for i=1:ENVCHANGES
            copyEnvState = copy(me.EnvState[:,envState]) # copy founder.EnvState without changing founder.EnvState
            copyEnvState[mutatedEnvs[mutatedEnvIndx[i]]] = 1 .- copyEnvState[mutatedEnvs[mutatedEnvIndx[i]]]
            mutatedEnvMatrix[:,i] = copyEnvState
        end
        return(mutatedEnvMatrix)
    end
end


# -------------------------------------------
# METHODS: Evolution Step
# -------------------------------------------
function update(mes::Vector{Individual{T}}) where T
# This function updates the state of a vector (Sw) of individuals
# ****One run of "update" function corresponds to one generation in evolution****

    # pmap runs in parallel if julia is invoked with multiple threads
    #pmap(develop, mes) # Run the "develop" function for the given input individual
    #*** Why running develop function for every generation (each iteration of update
    # function) because at end of while loop, will have replaced every individual
    # in the population and these individuals have undergo development during
    # while loop, so essentially developing individuals twice??? Except in the case
    # of the first generation development step - ML question 11/13/17
    oldinds = deepcopy(mes)

    newind = 1 # new individual who is stable and with
    # high enough fitness will replace the old individual

    while newind <= length(mes)
        #tempind = Individual()

        z = rand(1:N) # Pick a random individual from population
                    # with N individuals total to undergo asexual
                    # or sexual reproduction

        if SEXUALREPRO
            r = rand(1:N) # SEXUALREPRO is false for polycomb code and analysis --> not after 2018
            tempind = deepcopy(oldinds[z])
            reproduce(oldinds[z], oldinds[r], tempind)
        else
            # *** Why not just put tempind = copy(oldinds[z])?? -ML 4/17/18 --> fixed on 4/17/18
            tempind = deepcopy(oldinds[z]) #** is this the same as doing pop.individuals[z]? I don't think it is.... -ML 6/20/19 should input to update() be pop.individuals??? -ML 6/20/19
            #tempind.network = copy(oldinds[z].network)
            #tempind.initstate = copy(oldinds[z].initstate)
            #tempind.initstate2 = copy(oldinds[z].initstate2)
            #tempind.optstate = copy(oldinds[z].optstate)
            #tempind.optstate2 = copy(oldinds[z].optstate2)
            #tempind.develstate = copy(oldinds[z].develstate)
            # ** Why not saving .polycombstate as well?? -ML 4/17/18
            #tempind.polycombvec = copy(oldinds[z].polycombvec) # theta vector; representing polycomb susceptible genes, susceptible if = 1, not if = 0
            #tempind.EnvState = copy(oldinds[z].EnvState)
            #tempind.EnvIndInter = copy(oldinds[z].EnvIndInter)
        end

        mutate(tempind) # possibly mutate individuals gene network (W)

        if !SAURABVERSION & MUTATEENVINTERACTIONFLAG
            mutateEnvIndInter(tempind) # possibly mutate individuals environment gene network interaction matrix (W x E)
        end

        if POLYCOMBFLAG # True for polycomb code and analysis
            polymut(tempind)
        end

        iterateIndForOneEnvState(tempind, 1)

        if tempind.stable[1] # select for developmental stability for individual 1
            if DIFENVS > 1
                for j = 2:DIFENVS
                    iterateIndForOneEnvState(tempind, j)
                end

                if tempind.stable == trues(DIFENVS) # select for developmental stability for individual under all different environmental states
                    # If any is not stable, then don't keep
                    if MEASUREFIT
                        fitnesseval(tempind) # set to true in polycomb analysis
                        fitnessUnderEachEnvEval(tempind) # measure fitness under each environmental condition separately
                    else
                        tempind.fitness = 1.
                        tempind.fitnessUnderEachEnv = ones(DIFENVS)
                    end

                    if tempind.fitness >= rand() # Stabilizing selection step:
                            # higher fitness means more likely to
                            # have reproductive success because more likely to get picked
                            # out of population --> rand() produces numbers 0 to 1 with
                            # Uniform distribution, therefore the larger the number, the
                            # more likely it will be greater than the rand() number
                        mes[newind] = deepcopy(tempind)
                        newind += 1
                    end
                end
            end
        end
    end
end


# ------------------------------------------
# METHODS: Development Step
# ------------------------------------------

function develop(me::Individual)
# Run the normal developmental process for a single
# individual that assumes its environment state vector
# remains unchanged, and lets polycomb state vector
# evolve freely until critical time in development (TIMECRIT).
# Then sets polycomb state vector after TIMECRIT by
# completely activating polycomb response elements
# (setting polycomb state vector entry to 0) if the
# corresponding polycomb susceptible gene's expression
# is less than a certain threshold at TIMECRIT

    iterateind(me)

end


function iterateind(me::Individual)
# iterate individuals from their initial state to
# their developmental state and store the results
# in the Individual object
    currstate = copy(me.initstate) # currstate = Sw vector; initstate does NOT get changed
    stateupdate = zeros(Float64,G)# Define stateupdate & stateupdateenv before for loop because for
        # loop defines a new scope
    stateupdateenv = zeros(Float64,G)
    vectAD = zeros(Float64,MAXCONV) #vector to 'save' the actualDist
    if POLYCOMBFLAG # If zero then not running with no polycomb
        for i=1:TIMECRIT
            stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
            stateupdateenv = me.EnvIndInter*me.EnvState # (W x E)*Se: environment interaction matrix times environment state vector
            stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
            stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # converting resulting state vector into 0's or 1's --> if vector contains only -1's and 1's then will not alter vector
            currstate = copy(stateupdate)
            actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                             # and the current time point, Sw t+1 (stateupdate)
            vectAD[i] = actualDist
        end

        # run after timecrit
        polytemp = findall(x-> x < GENETHRESH, stateupdate) # find genes that are less than or equal to a gene threshold, GENETHRESH
        prcPresent = sum(me.polycombvec, dims = 2)
        prcPresentPos = findall(x -> x > 0, prcPresent)
        # to turn from Cartisean Index to list of indexes in vector:
        polys = map(x -> prcPresentPos[x][1], 1:length(prcPresentPos)) # Theta vector. Genes possibly susceptible to polycomb when
        # polycombvec = 1, 2, ..., PRCTOTAL: for example, gene has polycomb response element for either PRC1, PRC2 or both.
        # Thus, polycombstate entry needs to be set to 0 if this gene's expression is less than a gene threshold and
        # polycomb protein is active (assumed to be active until explicitly break PcG protein).
        polycombind = intersect(polytemp, polys) # Find index of genes that are susceptible to
        # polycomb elements (have a PRE & thus polycombvec >= 1) and less than or equal to a gene threshold after
        # critial time in development (TIMECRIT)
        #me.polycombstate = ones(Int64,G) # ** Why need this line of code during every step of evolution??? - ML 4/11/18 # set all polycomb susceptible genes as being expressed
        # (meaning polycomb response element (PcG protein) is inactive) bc set to 1
        me.polycombstate[polycombind] .= 0 # polycombstate entry is 0 if polycomb susceptible
        # gene is turned off by polycomb response element, i.e. the gene's corresponding
        # polycomb response element is active
    end

    for i=(1+TIMECRIT):MAXCONV
        stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
        stateupdateenv = me.EnvIndInter*me.EnvState # (W x E)*Se: environment interaction matrix times environment state vector
        stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
        stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # sigma s in eqs
        stateupdate = stateupdate .* me.polycombstate # turn off polycomb supressed genes
        # save vector of length TAU of stateupdate and currstate to input to slidingWindow function
        actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                         # and the current time point, Sw t+1 (stateupdate)
        # Sliding window - added 7/5/18 - SG
        # if number of iterations once testing for stability (so after TIMECRIT) is >= TAU,
        # then start calling sliding window function to test from (i-(TAU-1)) to i
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.))
            if JACOBIANFLAG
                jacobianMat = jacobian(me,stateupdate)
                eigenVals = eigvals(jacobianMat)
                largestEigVal = maximum(real(eigenVals))
                stableBool = (largestEigVal < 0.0)
            else # else run slidingWindow to test for stability
                stableBool = slidingWindow(i,vectAD)
            end
            if distBool && stableBool # returns true if
                me.stable = true
                me.develstate = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength = i #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable = false
                me.develstate = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength = MAXCONV #changed to copy() - SG 6/21/18
            end
        end

        currstate = copy(stateupdate)
    end
end


function iterateIndForOneEnvState(me::Individual, envState::Int64)
# iterate individuals from their initial state to
# their developmental state and store the results
# in the Individual object
    currstate = copy(me.initstate) # currstate = Sw vector; initstate does NOT get changed
    stateupdate = zeros(Float64,G)# Define stateupdate & stateupdateenv before for loop because for
        # loop defines a new scope
    stateupdateenv = zeros(Float64,G)
    vectAD = zeros(Float64,MAXCONV) #vector to 'save' the actualDist
    if POLYCOMBFLAG # If zero then not running with no polycomb
        for i=1:TIMECRIT
            stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
            stateupdateenv = me.EnvIndInter*me.EnvState[:,envState] # (W x E)*Se: environment interaction matrix times environment state vector
            stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
            stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # converting resulting state vector into 0's or 1's --> if vector contains only -1's and 1's then will not alter vector
            currstate = copy(stateupdate)
            actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                             # and the current time point, Sw t+1 (stateupdate)
            vectAD[i] = actualDist
        end

        # run after timecrit
        polytemp = findall(x-> x < GENETHRESH, stateupdate) # find genes that are less than or equal to a gene threshold, GENETHRESH
        prcPresent = sum(me.polycombvec, dims = 2)
        prcPresentPos = findall(x -> x > 0, prcPresent)
        # to turn from Cartisean Index to list of indexes in vector:
        polys = map(x -> prcPresentPos[x][1], 1:length(prcPresentPos)) # Theta vector. Genes possibly susceptible to polycomb when
        # polycombvec = 1, 2, ..., PRCTOTAL: for example, gene has polycomb response element for either PRC1, PRC2 or both.
        # Thus, polycombstate entry needs to be set to 0 if this gene's expression is less than a gene threshold and
        # polycomb protein is active (assumed to be active until explicitly break PcG protein).
        polycombind = intersect(polytemp, polys) # Find index of genes that are susceptible to
        # polycomb elements (have a PRE & thus polycombvec >= 1) and less than or equal to a gene threshold after
        # critial time in development (TIMECRIT)
        #me.polycombstate = ones(Int64,G) # ** Why need this line of code during every step of evolution??? - ML 4/11/18 # set all polycomb susceptible genes as being expressed
        # (meaning polycomb response element (PcG protein) is inactive) bc set to 1
        me.polycombstate[polycombind, envState] .= 0 # polycombstate entry is 0 if polycomb susceptible
        # gene is turned off by polycomb response element, i.e. the gene's corresponding
        # polycomb response element is active for given environmental state, envState (i.e. column of polycombstate)
    end

    for i=(1+TIMECRIT):MAXCONV
        stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
        stateupdateenv = me.EnvIndInter*me.EnvState[:,envState] # (W x E)*Se: environment interaction matrix times environment state vector
        stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
        stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # sigma s in eqs
        stateupdate = stateupdate .* me.polycombstate[:,envState] # turn off polycomb supressed genes
        # save vector of length TAU of stateupdate and currstate to input to slidingWindow function
        actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                         # and the current time point, Sw t+1 (stateupdate)
        # Sliding window - added 7/5/18 - SG
        # if number of iterations once testing for stability (so after TIMECRIT) is >= TAU,
        # then start calling sliding window function to test from (i-(TAU-1)) to i
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.))
            if JACOBIANFLAG
                jacobianMat = jacobian(me,stateupdate)
                eigenVals = eigvals(jacobianMat)
                largestEigVal = maximum(real(eigenVals))
                stableBool = (largestEigVal < 0.0)
            else # else run slidingWindow to test for stability
                stableBool = slidingWindow(i,vectAD)
            end
            if distBool && stableBool # returns true if
                me.stable[envState] = true
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = i #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable[envState] = false
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = MAXCONV #changed to copy() - SG 6/21/18
            end
        end

        currstate = copy(stateupdate)
    end
end


function reproduce(me::Individual, you::Individual, us::Individual) # Do not do sexual reproduction in polycomb stuff because individual cells behave more like asexual reproduction
# sexual reproduction via independent row segregation

    reproindxs = rand(0:1,G)

    for i in findall(x->x==1,reproindxs)
        us.network[i,:] = deepcopy(me.network[i,:])
    end

    for j in findall(x->x==0,reproindxs)
        us.network[j,:] = deepcopy(you.network[j,:])
    end

    return reproindxs
end


function mutate(me::Individual)
# mutate nonzero elements of an individuals
# network according to a rate parameter normalized
# by the size of the nonzero entries in an
# individual network
# done for each reproduction step for each generation

    # Find the non-zero entries as potential mutation sites
    nzindx = findall(x->x!=0,me.network)
    # Determine the connectivity of the matrix
    # by counting the number of non-zeros:
    cnum = length(nzindx) # note cnum=C*G^2, where c is the connectivity of network (W)
    mutflag = 0
    mutatedElement = 0
    for i=1:cnum # chance to mutate any of the non-zero elements in W matrix
    # with a total of 10% chance of a *single* mutation per individual
        if rand() < MUTRATE/cnum # about MUTRATE/cnum percent chance that
        # will mutate this particular non-zero element of W matrix
            me.network[nzindx[i]] = MUTMAG*randn() # Mutate a non-zero entry by a number chosen from a gaussian
            mutflag =+ 1
            mutatedElement = nzindx[i]
        end
    end
    return [mutflag, mutatedElement]
end


function mutateEnvIndInter(me::Individual)
# Function mutates nonzero elements of an individuals
# Environment Individual Interaction matrix (W x E)
# according to a rate parameter normalized
# by the size of the nonzero entries in an
# individual environment Interaction matrix.
# Function is ran for each reproduction step
# for each generation

    # Find the non-zero entries as potential mutation sites
    nzindx = findall(x->x!=0,me.EnvIndInter)

    # Determine the connectivity of the matrix
    # by counting the number of non-zeros
    cnum = length(nzindx)
    envmutflag=0
    mutatedElement=0
    for i=1:cnum
        # For each non-zero entry:
        # With probability Rate/C*G^2, note cnum=(C*(G^2))
        if  rand() < ENVMUTRATE/cnum
            # Mutate a non-zero entry by a number chosen from a gaussian
            me.EnvIndInter[nzindx[i]] = ENVMUTRATE*randn()
            envmutflag =+ 1
            mutatedElement = nzindx[i]
        end
    end
    return [envmutflag, mutatedElement]
end


function polymut(me::Individual)
# perform a one-site mutation of a given
# individual's polycomb vector
# to get population with individuals that can have
# different polycomb response elements inactivated
    mutmat = copy(me.polycombvec)
    for i = 1:G
        if rand() < POLYMUTRATE/G #10% of polycomb susceptible genes can be mutated and there are G polycomb susceptible genes
            prcToMutate = rand(1:PRCS) # select with PRC to mutate, so column to mutate because columns of polycombvec are PRCs when have more than one PRC
            mutmat[i,prcToMutate] = 1-mutmat[i,prcToMutate] # change zero entry to 1, or 1 entry to zero;
        end
    end
    me.polycombvec = mutmat
    return mutmat
end



#--------------------------------------------------------------
# Functions for breakage of polycomb and robustness testing
#--------------------------------------------------------------
function iterateind2(me::Individual, envState::Int64, envStateVecToTest::Vector{Float64})
# iterate individuals from their development state to
# their later state and store the results
# in the Individual object - version to be used with an already set polycomb state vector
# i.e. no iteration to TIMECRIT needed
# Testing for cell pliancy (so change EnvState but not polycomb state (so typically test
# using individual after went through development state)
    currstate = copy(me.initstate)
    stateupdate = zeros(Float64,G)# Define stateupdate & stateupdateenv before for loop because for
        # loop defines a new scope
    stateupdateenv = zeros(Float64,G)
    vectAD = zeros(Float64,MAXCONV)
    for i=1:MAXCONV
        stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
        stateupdateenv = me.EnvIndInter*envStateVecToTest # (W x E)*Se: environment interaction matrix times environment state vector
        stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
        stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # sigma s in eqs
        stateupdate = stateupdate .* me.polycombstate[:,envState] # turn off polycomb supressed genes
        # save vector of length TAU of stateupdate and currstate to input to slidingWindow function
        actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                         # and the current time point, Sw t+1 (stateupdate)
        # Sliding window - added 7/5/18 - SG
        # if number of iterations once testing for stability (so after TIMECRIT) is >= TAU,
        # then start calling sliding window function to test from (i-(TAU-1)) to i
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.))
            if JACOBIANFLAG
                jacobianMat = jacobian(me,stateupdate)
                eigenVals = eigvals(jacobianMat)
                largestEigVal = maximum(real(eigenVals))
                stableBool = (largestEigVal < 0.0)
            else # else run slidingWindow to test for stability
                stableBool = slidingWindow(i,vectAD)
            end
            if distBool && stableBool # returns true if
                me.stable[envState] = true
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = i #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable[envState] = false
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = MAXCONV #changed to copy() - SG 6/21/18
            end
        end
        currstate = copy(stateupdate)
    end
end

function iterateindForStability(me::Individual, envState::Int64, envStateVecToTest::Vector{Float64}, oldEnvState::Int64, indInOldEnv::Individual)
# iterate individuals from their development state to
# their later state and store the results
# in the Individual object - version to be used with an already set polycomb state vector
# i.e. no iteration to TIMECRIT needed --> so NOT going through development, just testing for stability
# Testing for cell pliancy (so change EnvState but not polycomb state (so typically test
# using individual after went through development state)
    currstate = copy(indInOldEnv.develstate[:,oldEnvState])
    stateupdate = zeros(Float64,G)# Define stateupdate & stateupdateenv before for loop because for
        # loop defines a new scope
    stateupdateenv = zeros(Float64,G)
    vectAD = zeros(Float64,MAXCONV)
    for i=1:MAXCONV
        stateupdate = me.network*currstate # W*Sw; gene interaction matrix times gene state vector
        stateupdateenv = me.EnvIndInter*envStateVecToTest # (W x E)*Se: environment interaction matrix times environment state vector
        stateupdate = stateupdate + stateupdateenv # W*Sw + (W x E)*Se
        stateupdate = 1.0./(1.0.+exp.(-SIGSTR*stateupdate)) # sigma s in eqs
        stateupdate = stateupdate .* me.polycombstate[:,oldEnvState] # turn off polycomb supressed genes that were turned off in old environment individual originally evolved in
        # save vector of length TAU of stateupdate and currstate to input to slidingWindow function
        actualDist = sum(abs.(currstate - stateupdate))/G # normalized Hamming distance between the previous time point state, Sw t (currstate),
                                                         # and the current time point, Sw t+1 (stateupdate)
        # Sliding window - added 7/5/18 - SG
        # if number of iterations once testing for stability (so after TIMECRIT) is >= TAU,
        # then start calling sliding window function to test from (i-(TAU-1)) to i
        vectAD[i] = actualDist
        if i >= TAU
            distBool = (vectAD[i] <= 10.0^(-6.))
            if JACOBIANFLAG
                jacobianMat = jacobian(me,stateupdate)
                eigenVals = eigvals(jacobianMat)
                largestEigVal = maximum(real(eigenVals))
                stableBool = (largestEigVal < 0.0)
            else # else run slidingWindow to test for stability
                stableBool = slidingWindow(i,vectAD)
            end
            if distBool && stableBool # returns true if
                me.stable[envState] = true
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = i #changed to copy() - SG 6/21/18
                break
            elseif i==MAXCONV
                me.stable[envState] = false
                me.develstate[:,envState] = copy(stateupdate) #changed to copy() - SG 6/21/18
                me.pathlength[envState] = MAXCONV #changed to copy() - SG 6/21/18
            end
        end

        currstate = copy(stateupdate)
    end
end

function breakPRC(me::Individual, prcToBreak::Array{Int64}, envState::Int64)
    # prcToBreak can be one or multiple Array element(s), where just one many PRC to break
    # (prcToBreak = [1]), OR can give multiple PRC's to break (prcToBreak = [1,3])
    polycombstateVec = copy(me.polycombstate[:,envState])
    polycombpos = findall(x->x == 0, polycombstateVec) # find polycomb response element on genes that are active (gene is supressed by polycomb), meaning polycombstate = 0
    polycombVec = copy(me.polycombvec)
    if length(polycombpos) > 0 # check that at least one polycomb response element is active before breaking one
        prcToBreakMat = polycombVec[:,prcToBreak] # pull out only polycombvec for prcToBreak
        prcToBreakPos = findall(prcToBreakMat.>0) # find all genes (positions) that are susceptible to give prcToBreak (polycombVec == 1 if susceptible)
        if length(prcToBreakPos) > 0 # If there is at least one gene susceptible to prcToBreak
            prcToBreakMat[prcToBreakPos] .= 0 # set these susceptible genes to no longer be susceptible, so set polycombvec to 0
            polycombVec[:,prcToBreak] = prcToBreakMat # update polycombvec for prcToBreak so that all susceptible genes are set to 0 so that no longer susceptible
            prcStillPresentVec = sum(polycombVec, dims = 2) # sum rows of polycombvec to find genes that are still susceptible to other active PRCs (prcStillPresentVec == 1 then gene susceptible to another active PRC)
            prcNotPresent = findall(x -> x == 0, prcStillPresentVec) # find position of genes that are NOT still susceptible to other active PRCs
            presToBreak = map(x -> prcNotPresent[x][1], 1:length(prcNotPresent)) # turn prcNotPresent from Cartisean Index to list of indexes in vector
            polycombstateVec[intersect(presToBreak, polycombpos)] .= 1 # set genes that were ONLY susceptible to prcToBreak as NOT being suppressed by that PRC (genes that were not suppressed don't need to have their polycombstate updated because they are already 1)
            me.polycombstate[:,envState] = polycombstateVec # update polycombstate for given envState with prcToBreak breakdown gene results
        else
            print("PRC(s) choosen to break not active")
        end
    else
        print("no PRC that is active")
    end
end

function breakpoly(me::Individual, envState::Int64)
# function takes an individual's existing polycomb state vector
# and abrogates an entry from 0 to 1 to make exactly one polycomb
# response element inactive (break it's polycomb such that
# polycomb can no longer supress polycomb susceptible gene)

# Original way to breakdown polycomb when only one universal PRC
    polycombstateVec = copy(me.polycombstate[:,envState])
    polycombpos = findall(x->x == 0, polycombstateVec) # find polycomb response element on genes that are active (gene is supressed by polycomb), meaning polycombstate = 0
    if length(polycombpos) > 0 # check that at least one polycomb response element is active before breaking one
        breakidx = rand(1:length(polycombpos))  # break exactly 1 polycomb response element, i.e. make the polycomb response element inactive by stating polycombstate = 1
        polycombstateVec[polycombpos[breakidx]] = 1
        me.polycombstate[:,envState] = polycombstateVec
    end
end

#--------------------------------------------------------------
# MEASUREMENTS FOR FITNESS, ROBUSTNESS, AND MODULARITY
#--------------------------------------------------------------
function fitnessevalWithoutAdjustment(me::Individual)
# measure fitness according to the distance
# between the developmental state determined
# by running iterateind(me) and the optimum
# state (founder development state) for the
# two individuals
    if me.stable == trues(DIFENVS)
        distanceForFitness = sum(map(x -> sum(abs.(me.develstate[:,x] - me.optstate[:,x])), 1:DIFENVS))/(G*DIFENVS) # ranges from 0 to infinity
        fitnessResult = exp.(-(distanceForFitness/SELSTR)) # me.fitness ranges from 1 to
        # zero, with 1 being the most fit (i.e. distance = 0)
        me.fitness = fitnessResult
    else
        me.fitness=0 # fitness set to zero when individual does not reach
        # developmental stability (steady state by MAXCONV iterations)
    end
end


function fitnesseval(me::Individual)
# measure fitness according to the distance
# between the developmental state determined
# by running iterateind(me) and the optimum
# state (founder development state) for the
# two individuals
    if me.stable == trues(DIFENVS)
        sumOfDistancesForFitness = sum(map(x -> sum(abs.(me.develstate[:,x] - me.optstate[:,x])), 1:DIFENVS))
        distanceForFitness = sumOfDistancesForFitness/(G*DIFENVS) # ranges from 0 to infinity
        fitnessResult = (exp.(-(distanceForFitness/SELSTR)) + maximum([0,(PEAK-(PEAK/MINFIT)*(sum(abs.(me.develstate[:,2] - me.optstate[:,2]))/(G)))]))/2 # me.fitness ranges from 1 to
        # zero, with 1 being the most fit (i.e. distance = 0)
        me.fitness = fitnessResult
    else
        me.fitness=0 # fitness set to zero when individual does not reach
        # developmental stability (steady state by MAXCONV iterations)
    end
end


function fitnessUnderEachEnvEval(me::Individual)
    if me.stable == trues(DIFENVS)
        for i = 1:DIFENVS
            distanceForFitness = sum(abs.(me.develstate[:,i] - me.optstate[:,i]))/G # ranges from 0 to infinity
            fitnessResult = exp.(-(distanceForFitness/SELSTR)) # me.fitness ranges from 1 to
            # zero, with 1 being the most fit (i.e. distance = 0)
            me.fitnessUnderEachEnv[i] = fitnessResult
        end
    else
        me.fitnessUnderEachEnv=zeros(Float64,DIFENVS) # fitness set to zero when individual does not reach
        # developmental stability (steady state by MAXCONV iterations)
    end
end

function founderGeneticRobust(me::Individual)
# measure sensitivity to mutations in gene network (W matrix)
# of the developmental state (phenotypic output vector S-hat)
    local tempdiff
    stablePerturbations = zeros(Int64,DIFENVS) # for unit testing purposes
    for j = 1:DIFENVS
        dist = 0.
        for i=1:10000
            perturbed = deepcopy(me)
            perturbed.network = onemut(me) # mutate one non-zero element in the W matrix (gene interaction network)
            iterateIndForOneEnvState(perturbed, j) # test for stability
            if perturbed.stable[j] # for unit testing purposes
                stablePerturbations[j] = stablePerturbations[j] + 1 # for unit testing purposes
            end
            tempdiff = me.develstate[:,j] - perturbed.develstate[:,j] # should we only measure robustness for ones that are stable?? *ML 5//29/19
            dist += sum(abs.(tempdiff))/(G)
        end
        me.robustness[j] = dist/10000 # resulting average genetic robustness for
        # number of times introduce one mutation (ROBIT)
    end
    return stablePerturbations # for unit testing purposes
end

function geneticRobustness(me::Individual, founderGeneticRobustness::Vector{Float64})
# measure sensitivity to mutations in gene network (W matrix)
# of the developmental state (phenotypic output vector S-hat)
# you is founder individual
    local tempdiff
    stablePerturbations = zeros(Int64,DIFENVS) # for unit testing purposes
    for j = 1:DIFENVS
        dist = 0.
        for i=1:ROBIT # ROBIT is # of different mutations testing
            perturbed = deepcopy(me)
            perturbed.network = onemut(me) # mutate one non-zero element in the W matrix (gene interaction network)
            iterateIndForOneEnvState(perturbed, j) # test for stability
            if perturbed.stable[j] # for unit testing purposes
                stablePerturbations[j] = stablePerturbations[j] + 1 # for unit testing purposes
            end
            tempdiff = me.develstate[:,j] - perturbed.develstate[:,j] #**Should this be
            # me.develstate or founder.develstate?? Was initially me.develstate.
            # Aviv thinks it should be founder.develstate to match Saurab model -ML 05/22/18
            # If selection strength is very high, then me.develstate will be equal to founder.develstate
            dist += sum(abs.(tempdiff))/(G)
        end
        me.robustness[j] = founderGeneticRobustness[j] - (dist/ROBIT)
    end
    return stablePerturbations # for unit testing purposes
end

function polycombGeneticRobustness(me::Individual, founderGeneticRobustness::Vector{Float64})
# measure sensitivity to mutations in gene network (W matrix)
# of the developmental state (phenotypic output vector S-hat) when already
# have set polycomb state vector. Difference with genetic robustness function
# is that do not set robustness value to the individual testing, just return the
# robustness value because this function is for testing purposes with polycomb.
    local tempdiff
    for j = 1:DIFENVS
        dist = 0.
        for i=1:ROBIT
            perturbed = deepcopy(me)
            perturbed.initstate = copy(perturbed.develstate[:,j]) # set the initial state
                # of the individual to that of the development state (Sw) at the end
                # of evolution
            perturbed.network = onemut(me)
            iterateIndForOneEnvState(perturbed, j) # use this function since already have set polycomb vector
            tempdiff = me.develstate[:,j] - perturbed.develstate[:,j]
            dist += sum(abs.(tempdiff))/(G)
        end
        robustness[j] = founderGeneticRobustness[j] - (dist/ROBIT)
    end
    return robustness
end

function founderEnvRobust(me::Individual)
# measure sensitivity of Founder to changes in environmental conditions
# present of the developmental state (phenotypic output vector S-hat w)
    local tempdiff
    stablePerturbations = zeros(Int64,DIFENVS) # for unit testing purposes
    for j = 1:DIFENVS
        dist = 0.0
        if SAURABVERSION
            mutatedInitStates = collect(combinations(collect(1:length(me.initstate)), convert(Int64, ALTERENV*G)))
            mutatedInitStateMatrix = zeros(Float64, G, length(mutatedInitStates))
            for i=1:length(mutatedInitStates)
                copyInitState = copy(me.initstate) # copy founder.initstate without changing founder.initstate or me.initstate
                copyInitState[mutatedInitStates[i]] = 1 .- copyInitState[mutatedInitStates[i]]
                mutatedInitStateMatrix[:,i] = copyInitState
            end
            for i=1:length(mutatedInitStates)
                perturbed = deepcopy(me)
                perturbed.initstate = copy(mutatedInitStateMatrix[:,i])
                iterateIndForOneEnvState(perturbed, j) # use this function to test if stability is reached
                if perturbed.stable[j] # for unit testing purposes
                    stablePerturbations[j] = stablePerturbations[j] + 1 # for unit testing purposes
                end
                tempdiff = me.develstate[:,j] - perturbed.develstate[:,j]
                dist += sum(abs.(tempdiff))/(G)
            end
            me.envRobustness[j] = dist/length(mutatedInitStates)
        else
            mutatedEnvs = collect(combinations(collect(1:length(me.EnvState[:,j])), convert(Int64, ALTERENV*ENVS))) # all combinations of
                # changing only 20% (ALTERENV) of the elements of Environment state vector (ENVS)
            mutatedEnvMatrix = zeros(Float64, ENVS, length(mutatedEnvs))
            for i=1:length(mutatedEnvs)
                copyEnvState = copy(me.EnvState[:,j]) # copy founder.EnvState without changing founder.EnvState
                copyEnvState[mutatedEnvs[i]] = 1 .- copyEnvState[mutatedEnvs[i]]
                mutatedEnvMatrix[:,i] = copyEnvState
            end
            for i=1:length(mutatedEnvs)
                perturbed = deepcopy(me)
                perturbed.EnvState[:,j] = copy(mutatedEnvMatrix[:,i])
                iterateIndForOneEnvState(perturbed, j) # use this function to test if stability is reached
                if perturbed.stable[j] # for unit testing purposes
                    stablePerturbations[j] = stablePerturbations[j] + 1 # for unit testing purposes
                end
                tempdiff = me.develstate[:,j] - perturbed.develstate[:,j]
                dist += sum(abs.(tempdiff))/(G)
            end
            me.envRobustness[j] = dist/length(mutatedEnvs)
        end
    end
    return stablePerturbations # for unit testing purposes
end

function environmentalRobustness(me::Individual, mutatedEnvMatrix::Array{Array{Int64,2},1}, founderEnvironmentalRobustness::Vector{Float64})
# measure sensitivity of each individual in population to changes in
# environmental conditions present of the developmental state (phenotypic
# output vector S-hat w). you is founder individual.
    local tempdiff
    stablePerturbations = zeros(Int64,DIFENVS) # for unit testing purposes
    for j = 1:DIFENVS
        dist = 0.0
        if SAURABVERSION
            for i=1:ENVCHANGES
                perturbed = deepcopy(me)
                perturbed.initstate = copy(mutatedEnvMatrix[j][:,i])
                iterateIndForOneEnvState(perturbed, j) # use this function to test if stability is reached
                if perturbed.stable[j] # for unit testing purposes
                    stablePerturbations[j] = stablePerturbations[j] + 1 # for unit testing purposes
                end
                tempdiff = me.develstate[:,j] - perturbed.develstate[:,j] #**Should this be
                # me.develstate or founder.develstate?? Was initially me.develstate.
                # Aviv thinks it should be founder.develstate to match Saurab model -ML 05/22/18
                dist += sum(abs.(tempdiff))/(G)
            end
        else
            for i=1:ENVCHANGES
                perturbed = deepcopy(me)
                perturbed.EnvState[:,j] = copy(mutatedEnvMatrix[j][:,i])
                iterateIndForOneEnvState(perturbed, j) # use this function to test if stability is reached
                if perturbed.stable[j] # for unit testing purposes
                    stablePerturbations[j] = stablePerturbations[j] + 1 # for unit testing purposes
                end
                tempdiff = me.develstate[:,j] - perturbed.develstate[:,j]
                dist += sum(abs.(tempdiff))/(G)
            end
        end
        # Calculate environmental robustness by subtracting individual's
        # sensitivity to environmental changes from the Founder's sensitivity
        # to environmental changes:
        me.envRobustness[j] = founderEnvironmentalRobustness[j] - (dist/ENVCHANGES)
    end
    return stablePerturbations # for unit testing purposes
end

function onemut(me::Individual)
# Function for robustness testing:
# perform a one-site mutation of a given
# individual's network
# This function does NOT alter me.network,
# and instead returns mutated network as "mutmat"
    nzindx = findall(x->x!=0,me.network) # find all non-zero elements in gene network (W) matrix
    cnum = length(nzindx) # number of non-zero elements in gene network
    i=rand(1:cnum) # randomly pick one element in network to mutate
    mutmat = copy(me.network)
    mutmat[nzindx[i]] = MUTMAG*randn()
    return mutmat
end


### BELOW NOT UPDATED FOR EVOLVING IN MULTIPLE ENVIRONMENTS
function environmentalRobustnessWithoutSelection(me::Individual, mutatedEnvMatrix::Matrix, founderEnvironmentalRobustness::Float64)
# measure sensitivity of each individual in population to changes in
# environmental conditions present of the developmental state (phenotypic
# output vector S-hat w). Compare to environmental robustness when
# do not perturb environment.
    local tempdiff
    stablePerturbations = 0 # for unit testing purposes
    dist = 0.0
    if SAURABVERSION
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.initstate = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            if perturbed.stable # for unit testing purposes
                stablePerturbations = stablePerturbations + 1 # for unit testing purposes
            end
            tempdiff = me.develstate - perturbed.develstate #**Should this be
            # me.develstate or founder.develstate?? Was initially me.develstate.
            # Aviv thinks it should be founder.develstate to match Saurab model -ML 05/22/18
            dist += sum(abs.(tempdiff))/(G)
        end
    else
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.EnvState = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            if perturbed.stable # for unit testing purposes
                stablePerturbations = stablePerturbations + 1 # for unit testing purposes
            end
            tempdiff = me.develstate - perturbed.develstate
            dist += sum(abs.(tempdiff))/(G)
        end
    end
    # Calculate environmental robustness by subtracting individual's
    # sensitivity to environmental changes from the Founder's sensitivity
    # to environmental changes:
    me.envRobustness = me.envRobustness - (dist/ENVCHANGES)
    return stablePerturbations # for unit testing purposes
end


function environmentalRobustnessSaurabhCcodeVersion(me::Individual, mutatedEnvMatrix::Matrix, founderEnvironmentalRobustness::Float64, founder::Individual)
# measure sensitivity of each individual in population to changes in
# environmental conditions present of the developmental state (phenotypic
# output vector S-hat w). you is founder individual.
    local tempdiff
    dist = 0.0
    if SAURABVERSION
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.initstate = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            tempdiff = founder.develstate - perturbed.develstate # founder.develstate to match Saurab model
            dist += sum(abs.(tempdiff))/(G)
        end
    else
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.EnvState = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            tempdiff = founder.develstate - perturbed.develstate
            dist += sum(abs.(tempdiff))/(G)
        end
    end
    # Calculate environmental robustness by subtracting individual's
    # sensitivity to environmental changes from the Founder's sensitivity
    # to environmental changes:
    me.envRobustness = founderEnvironmentalRobustness - (dist/ENVCHANGES)
end


function environmentalSensitivity(me::Individual, mutatedEnvMatrix::Matrix)
# measure sensitivity of each individual in population to changes in
# environmental conditions present of the developmental state (phenotypic
# output vector S-hat w). you is founder individual.
    local tempdiff
    dist = 0.0
    if SAURABVERSION
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.initstate = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            tempdiff = me.develstate - perturbed.develstate #**Should this be
            # me.develstate or founder.develstate?? Was initially me.develstate.
            # Aviv thinks it should be founder.develstate to match Saurab model -ML 05/22/18
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
    # Calculate environmental robustness by subtracting individual's
    # sensitivity to environmental changes from the Founder's sensitivity
    # to environmental changes:
    me.envRobustness = dist/ENVCHANGES
end


function environmentalSensitivitySaurabhCcodeVersion(me::Individual, mutatedEnvMatrix::Matrix, founder::Individual)
# measure sensitivity of each individual in population to changes in
# environmental conditions present of the developmental state (phenotypic
# output vector S-hat w). you is founder individual.
    local tempdiff
    dist = 0.0
    if SAURABVERSION
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.initstate = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            tempdiff = founder.develstate - perturbed.develstate #**Should this be
            # me.develstate or founder.develstate?? Was initially me.develstate.
            # Aviv thinks it should be founder.develstate to match Saurab model -ML 05/22/18
            dist += sum(abs.(tempdiff))/(G)
        end
    else
        for i=1:ENVCHANGES
            perturbed = deepcopy(me)
            perturbed.EnvState = copy(mutatedEnvMatrix[:,i])
            iterateind(perturbed) # use this function to test if stability is reached
            tempdiff = founder.develstate - perturbed.develstate
            dist += sum(abs.(tempdiff))/(G)
        end
    end
    # Calculate environmental robustness by subtracting individual's
    # sensitivity to environmental changes from the Founder's sensitivity
    # to environmental changes:
    me.envRobustness = dist/ENVCHANGES
end
