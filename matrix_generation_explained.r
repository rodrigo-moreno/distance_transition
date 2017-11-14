Implicit.InterAttractor.Simulation <- function(Network, P.error, Nreps) {
  # Simulates binomial mutations in all space
  #Get State Space and x(t+1) = f(x) for each state
  AttrsLandscape <- Get.Attractors.Landscape(Network)		### ALL POSSIBLE STATES + THEIR BASIN OF ATTRACTION
  Ngenes <- ncol(AttrsLandscape[[1]])
  Nattractors <- length(names(table(AttrsLandscape[[2]])))
  StateSpace <- AttrsLandscape[[1]]				### AGAIN, ALL POSSIBLE STATES
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
  ### SO, THE PREVIOUS LINE GIVES A MATRIX WITH A TON OF ELEMENTS, EITHER "1" OR "0", WITH A PROBABILITY OF BEING 1 OF P.error.
  
  #Create concatenated vector of X(t+1)
  NextSs <- as.numeric(apply(Next.T.StateSpace,1, function(i) rep(i, Nreps)))

  #Simulate "errors" in X(t+1)
  Mutind <- which(MutMatrix==1)				### WHICH ELEMENTS OF MutMatrix ARE 1'S?
  ZeroInd <- Mutind[which(NextSs[Mutind]==1)]		### THE NEXT FEW LINES CHANGE THE VALUE OF THE NEXT STATE TO THE OPPOSITE (FROM "0" TO "1" AND VICE VERSA) ON
  UnoInd <- Mutind[which(NextSs[Mutind]==0)]		### THE POSITIONS WHERE THERE WAS A MUTATION
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
