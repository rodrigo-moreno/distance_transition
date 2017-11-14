setwd("C:/Users/Rodrigo/Documents/Maestria/toygraph")
library(BoolNet)

### GET DISTANCE BETWEEN ATTRACTORS
g = loadNetwork("toygraph.txt")
att = getAttractors(g)
att_num = max(att$stateInfo$attractorAssignment)
atts = list()
for (i in c(1:att_num)){
  a = unlist(getAttractorSequence(att,i))
  print(a)
  print(typeof(a))
  atts[[i]] = as.vector(a)
}

Dist = matrix(0, att_num, att_num)
names = c(1:att_num)
rownames(Dist) = names
colnames(Dist) = names
### THE DISTANCE BETWEEN TWO POINTS IS THE SQRT OF THE SUM OF SQUARES OF THE DIFFERENCE IN EACH COORDINATE
for (i in c(1:att_num)){
  for (j in c((i+1):att_num)){
    if (j == 5){
      break
    }
    att1 = unlist(atts[i])
    att2 = unlist(atts[j])
    diff = att1 - att2
    diff = diff**2
    d = sqrt(sum(diff))
    Dist[j,i] = d
  }
}
Dist = as.dist(Dist)



### GET TRANSITION PROBABILITY FROM ONE ATTRACTOR TO THE OTHER ###
Get.Attractors.Landscape <- function(Network) {
  ### RETURNS A LIST WITH TWO ELEMENTS:
  ### 1. ALL INITITATION STATES FOR THE NETWORK
  ### 2. THE ATTRACTOR BASIN TO WHICH SAID INITIAL STATE BELONGS TO
  attrs <- getAttractors(Network)
  TransTable <- getTransitionTable(attrs)

  Inits <- TransTable[[1]]
  for(i in 2:length(Network$genes)) Inits <- cbind(Inits, TransTable[[i]])

  AttractsClusterVector <- TransTable$attractorAssignment

  return(list(Inits, AttractsClusterVector))
}

Implicit.InterAttractor.Simulation <- function(Network, P.error, Nreps) {
  # Simulates binomial mutations in all space
  #Get State Space and x(t+1) = f(x) for each state
  AttrsLandscape <- Get.Attractors.Landscape(Network)
  Ngenes <- ncol(AttrsLandscape[[1]])
  Nattractors <- length(names(table(AttrsLandscape[[2]])))
  StateSpace <- AttrsLandscape[[1]]
  Next.T.StateSpace <- StateSpace*0
  Next.T.StateSpace <- t(sapply(1:nrow(Next.T.StateSpace), function(i) Next.T.StateSpace[i,] <- stateTransition(Network, StateSpace[i,]))) ### THIS GETS THE NEXT STATE FOR EVERY STATE
  colnames(StateSpace) <- colnames(Next.T.StateSpace)
  Char.State.Space <- apply(StateSpace, 1, function(i) paste(i, collapse=""))

  #Create Attractors Transition Probability Matrix
  T.Prob.Mat <- matrix(0, Nattractors, Nattractors)
  rownames(T.Prob.Mat) <- 1:Nattractors
  colnames(T.Prob.Mat) <- 1:Nattractors
  AttrsInd <- as.numeric(colnames(T.Prob.Mat))

  #Create Muation Indicator Vector
  MutMatrix <- rbinom(Nreps*Ngenes*nrow(StateSpace), 1, P.error)  ### rbinom(n, s, p) GIVES A RANDOM DISTRIBUTION WITH n OBSERVATIONS, s MAX VALUE AND p PROBABILITY.

  #Create concatenated vector of X(t+1)
  NextSs <- as.numeric(apply(Next.T.StateSpace,1, function(i) rep(i, Nreps)))

  #Simulate "errors" in X(t+1)
  Mutind <- which(MutMatrix==1)
  ZeroInd <- Mutind[which(NextSs[Mutind]==1)]
  UnoInd <- Mutind[which(NextSs[Mutind]==0)]
  NextSs[ZeroInd] <- 0
  NextSs[UnoInd] <- 1

  #Split and count states in X(t+1). Match them with basins.
  NextSs <- apply(matrix(NextSs, Nreps*nrow(StateSpace), Ngenes, byrow=TRUE), 1, function(i) paste(i, collapse=""))
  NextSs.LL <- lapply(split(NextSs, rep(1:nrow(StateSpace), each=Nreps)), table)
  Basins.L <-  lapply(NextSs.LL, function(i) AttrsLandscape[[2]][match(names(i), Char.State.Space)])

  #Update Attractors Transition Probability Matrix
  for(j in 1:nrow(StateSpace)) T.Prob.Mat[AttrsLandscape[[2]][j], ] <- T.Prob.Mat[AttrsLandscape[[2]][j], ] + sapply(1:length(AttrsInd), function(i) sum(NextSs.LL[[j]][which(AttrsInd[i]==Basins.L[[j]])]))

  #Normalize Attractors Transition Probability Matrix
  T.Prob.Mat <- t(sapply(1:nrow(T.Prob.Mat), function(i) T.Prob.Mat[i,]/sum(T.Prob.Mat[i,])))
  return(T.Prob.Mat)
}


Pi = Implicit.InterAttractor.Simulation(g, 0.2, 100)  %%% I STILL HAVE TO SEE HOW THIS CORRELATES TO Dist AND BASIN SIZE WITH DIFFERENT VALUES OF P.error


### ONCE I OBTAIN Pi FOR DIFFERENT VALUES OF P.error, I CAN CORRECT THE 
### TRANSITION PROBABILITIES ACCORDING TO THE ATTRACTOR SIZE
### THEN I CAN COMPARE Pi WITH THE DISTANCES
corr = matrix(0, att_num, att_num)
for (i in c(1:att_num)){
  basin_size = att$attractors[[i]]$basinSize
  corr[i, 1:4] = basin_size
}
Pi_c = Pi * corr
