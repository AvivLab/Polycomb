# Testing adding two environments to evolve with to try to evolve such that reach output that is 40% different.
# By Maryl Lambros, August 2019

## Test for finding envState2 that is stable for an individual but 70% different than envState
## assuming that initstate1 are the same for both environmental state cases
pop = genpop()
popEnvState2Test = deepcopy(pop)
popEnvState2Test.individuals[1].stable = false
counts = 1
while popEnvState2Test.individuals[1].stable == false && counts < MAXENVINPUT
    EnvState2 = copy(pop.individuals[1].EnvState)
    EnvState2ToChange = randperm(ENVS)[1:convert(Int64,PERCENTDIFENVS*ENVS)]
    EnvState2[EnvState2ToChange] = 1 .- pop.individuals[1].EnvState[EnvState2ToChange]
    popEnvState2Test.individuals[1].EnvState = copy(EnvState2)
    iterateind(popEnvState2Test.individuals[1])
    #popEnvState2Test.individuals[1].stable
    global counts += 1
end
print("\n",counts)
distance = sum(abs.(pop.individuals[1].develstate - popEnvState2Test.individuals[1].develstate))/G

## Test for finding Sw2opt (optstate2) that is at least 40% different than Sw1opt (opstate)
distanceOptStates = 0.0
countNum = 1
while distanceOptStates < PERCENTDIFOPTSTATES && countNum < MAXOPTSTATES
    global optState2 = copy(pop.individuals[1].optstate)
    optState2ToChange = randperm(G)[1:convert(Int64,PERCENTDIFGENES*G)]
    optState2[optState2ToChange] = abs.(rand() .- pop.individuals[1].optstate[optState2ToChange])
    global distanceOptStates = sum(abs.(pop.individuals[1].optstate - optState2))/G
    global countNum += 1
end
print("\n",distanceOptStates)
print("\n",countNum)
popEnvState2Test.individuals[1].optstate = optState2

## Test for calculating fitness:
distanceForFitness = (sum(abs.(pop.individuals[1].develstate - pop.individuals[1].optstate)) + sum(abs.(popEnvState2Test.individuals[1].develstate - popEnvState2Test.individuals[1].optstate)))/G
fitnessTest = exp.(-(distanceForFitness/SELSTR))
