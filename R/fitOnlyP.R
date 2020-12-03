bayes_model <- "

var mixmu[Nn], mixtau[Nn], mixpi[Nn];

data
{
## Turning the input standard errors into precisions for the estimated log geometric mean product
## concentrations
  for (i in 1:Mn)
    {
      tau.se[i] <- 1/(se[i]*se[i])
    }
}

model
{
  ## +++++++++++++++++++++++++++++++++++++++++++++
  ## The heart of the model:
  ## Model the parent exposures
  ## +++++++++++++++++++++++++++++++++++++++++++++
  for (i in 1:N) {
    ## log of product of production volume and unit exposure
    lP[i] ~ dnorm(lPmu, tau.V)
    P[i]<- exp(lP[i])
  }
  ## hyperparameter for lP
  lPmu ~ dnorm(0, 0.001)

  ## code so that prior for sd(log(P[i])) is half-Cauchy(25)
  sd.dum ~ dt(0,1/625,1)
  sd.V <- abs(sd.dum)
  tau.V <- 1/pow(sd.V,2)
  ## ++++++++++++++++++++++++++++++++++++++++++++++
  ## Link the unobserved parent exposure to the observed
  ## (but censored) metabolite (w/some parent) urinary
  ## values.
  ##
  ## Model metabolite concentrations as (usually) unknown fractions
  ## of the parent.
  ##
  ## Phi[i,j] say what fraction of parent[i] goes to metabolite[j]
  ## observed values are lognormally distributed around U, with known
  ## precision, but censored.  Some of the elements of Phi are known to be
  ## 0 or 1 (thanks to Jim Rabinowitz for some of the 1's).
  lU <- log(t(Phi) %*% P)
  ## tau.se <- 1/(se * se)
  for (j in 1:Mn) {
    ly[j] ~ dnorm(lU[j],tau.se[j])
  }
  ## for j > Mn, what we have is the number of observations > LOD, and the total sample size.
  ## so let Pralod[j]  <-  1 - pnorm(lu[j], sd[j - Mn])
  ## Then ly[j]  <-  number above lod, and is binomial with mean Pralod[j] and n SS[j]
  for (j in (Mn+1):M) {
    Pralod[j - Mn] <- 1 - pnorm(lod[j - Mn], lU[j], tau[j - Mn])
    ly[j] ~ dbin(Pralod[j - Mn], SS[j - Mn])
  }
  ## The population lsds look to be mixtures of normal distributions for the data
  ## where moments seem well-estimated.
  ## Estimate mixmu, mixtau, mixpi externally, and input as data.
  for (j in 1:(M - Mn)) {
    lsd[j] ~ dnormmix(mixmu, mixtau, mixpi)
    tau[j]  <- exp(-2*lsd[j])
  }

  for (i in 1:NBranches)
  {
    phi[Bstart[i]:Bstop[i]] ~ ddirch(Alpha[Bstart[i]:Bstop[i]])
  }

  for (i in 1:Ndelta) {
    Phi[indx[i,1],indx[i,2]] <- phi[i]
  }

}
"


#' fitOnlyP
#'
#' Runs Bayesian inference for a provided subpopulation
#'
#' @param SUBPOP The subpopulation name for which we run the inference
#' @param Measured Output from the examine_error function (saved as
#'                 NewMeasured), with the last column being PosteriorShape
#' @param mapping  Output from the readNHANES function, saved in MCMCdata file
#' @param pred.data Output from the readNHANES function, saved in MCMCdata file
#' @param quick Default is FALSE.  If TRUE, the calculations will be run
#' with less iterations, resulting in a shorter overall runtime, however,
#' model convergence may not be reached.  Check convergence using the fromMolar function.
#' @param cores Number of cores to use during computation
#'
#' @import doMC
#' @import foreach
#' @import rjags
#' @import random
#' @import MASS
#' @import nor1mix
#' @importFrom coda varnames
#'
#' @return out.samps3R
#'         out.coda3R
#'
#' @export
#'
#'
fitOnlyP <- function(SUBPOP, Measured, mapping, pred.data, quick = FALSE, cores = 3) {

  registerDoMC(cores=cores)
  load.module("lecuyer")
  load.module("mix")
  doONEchain <- TRUE

  ## JAGS model is in OnlyP.jag
  ## It needs the following data:
  ##   y : observed metabolite rates (NA for below LOD) (1:M),
  ##   N = number of parents,
  ##   N.nontriv = number of parents with > 1 metabolite (sorted so they come first
  ##     in R.x)
  ##   N.mets[1:N.nontriv] = number of metabolites for parent i
  ##   indx[N,M]: indx[i,1:Nmets[i]] gives the indices in U for the metabolites
  ##     i contributes.  Corresponding values of Phi[i,] are non-zero.
  ##     If i > N.nontriv, indx[i,1] gives the entry in Phi[i,] that is 1
  ##     For j > N.nontriv, indx[i,j] gives column indices in Phi[i,] for 0's.
  ##   M = number of metabolites,
  ##   tau[1:M] : precision for observed metabolite rates.
  ##   lod = vector of limits of detection, tau = vector of precisions for observed
  ##   metabolite rates.
  ##   dstart[1:N.nontriv], dstop[1:N.nontriv]: indices into delta and phi
  ##     such that, e.g., delta[dstart[i]:dstop[i]] are the deltas for row i of Phi.
  ##     dstart[1] == 1, dstop[N.nontriv] == Ndelta.
  ##
  ## The following variables need initializing
  ##   P[i], i in 1:N
  ##   X1
  ##   X2
  ##   delta[1:Ndelta]: the non-trivial deltas corresponding to the non-trivial phis

  ## Data are based on John W.'s RData file DataTables-040212.RData, and successors.

  Measured <- subset(Measured, subpop == SUBPOP)
  M <- nrow(Measured)
  Mn <- sum(Measured$PosteriorShape == "normal")

  ## Now tailor the mapping information and pred.data to
  ## match the actual measurements we have.  Easiest to just
  ## recreate Phi and pred.data

  ## Drop entries for metabolites we won't be following: Drop CAS.1 that are NOT in Measured$CAS
  indx <- which(!(mapping$CAS.1 %in% unique(Measured$CAS)))
  if (length(indx) > 0)
    mapping <- mapping[-indx,]

  ## Sort this on CAS and Branch
  indx <- order(mapping$CAS, mapping$Branch)
  mapping <- mapping[indx,]

  ## Drop entries in pred.data for which we have lost Measured info.
  todrop <- setdiff(pred.data$CAS, mapping$CAS)
  DROP <- pred.data$CAS %in% todrop
  pred.data <- pred.data[!DROP,]

  ## Now, set up Phi to reflect this information
  ## -----------------------------------------------------------------------------
  ##                          Xref
  ##
  ## Build a matrix with rows for parents and columns for metabolites
  ## that is, rows for the records in pred.data, and columns for the
  ## records in Measured.  use CAS numbers to do the linking.
  ## The entries will be NA for values that will be estimated,
  ## 1 for entries representing 1:1 molar relationships between parent and metabolite,
  ## and 0 everywhere else (the majority of entries).
  ## Set all the entries to 0, initially.
  ## Then, set the cells corresponding to entries in NHANES.parent.metabolite.mapping
  ## to NA.  Some of these will be reset to 1,
  ## and the directly observed parents.  These are the CAS numbers in pred.data that
  ## are also in NHANES.
  Xref <- matrix(0, nrow=nrow(pred.data), ncol=length(unique(Measured$CAS)),
                 dimnames=list(pred.data$CAS, unique(Measured$CAS)))
  indx <- cbind(match(mapping$CAS,rownames(Xref)),
                match(mapping$CAS.1,colnames(Xref)))
  Xref[indx] <- NA

  Brnch <- by(mapping, mapping$CAS,
              FUN=function(z)
              {
                by(z, z$Branch, FUN=function(w) w$CAS.1, simplify=FALSE)
              }, simplify=FALSE
  )
  ## Brnch has an element for each row in Xref.  Each element has an element for each branch.
  ## Each branch with only one element corresponds to a 1 in Xref.  Each branch with multiple elements
  ## corresponds to that many rows in F with 1's in the columns corresponding to the multiple elements.

  for(pname in names(Brnch))
  {
    zz <- Brnch[[pname]]
    for (brn in names(zz))
    {
      if (length(zz[[brn]]) == 1)
      {
        Xref[ pname,zz[[brn]] ] <- 1
        Brnch[[pname]][brn] <- NULL
      }
    }
    if (length(Brnch[[pname]]) == 0) Brnch[pname] <- NULL
  }

  ## The remaining NAs in Xref correspond to deltas that need to be estimated.
  ## Let Ndeltas be the number of deltas:
  Ndelta <- sum(is.na(Xref))

  ## index[i,1], index[i,2] gives the row and column indices into Xref of the
  ## ith delta (or phi).  Make sure the rows and columns of Xref don't change
  ## order after index is created, and that the order of rows in Measured and
  ## pred.data match the order in Xref.

  ## Resort Measured so that the normal posteriors come first, then
  ## the uniform ones.  Break ties by sorting on CAS
  indx <- order(Measured$PosteriorShape, Measured$CAS)
  Measured <- Measured[indx,]
  ## Now, reorder the columns of Phi
  Xref <- Xref[,Measured$CAS]
  pred.data <- pred.data[rownames(Xref),]

  ## Now, loop through Brnch, fill index, and set up F
  index <- matrix(0,nrow=Ndelta, ncol=2)
  iF <- 1
  for(pname in names(Brnch))
  {
    i <- match(pname, rownames(Xref))           # location of this parent in Xref
    zz <- Brnch[[pname]]                        # metabolites for this parent
    for (brn in names(zz))                      # looping through metabolites
    {
      J <- match(zz[[brn]], colnames(Xref))     # columns in Xref
      for (k in 1:length(J))                    # loop through columns indices of Xref
      {
        K <- k + iF - 1                         # K is a counter for all of these non 1-to-1 relations
        index[K,1] <- i                         # row index of Xref for parent
        index[K,2] <- J[k]                      # column indices of Xref for metabolites
      }
      iF <- iF + length(J)
    }
  }
  nms <- paste(index[,1],index[,2],sep=":")
  NBranches <- sum(sapply(Brnch, length))
  Bstart <- Bstop <- numeric(NBranches)
  iB <- 0
  iF <- 1
  for(pname in names(Brnch))                    # loop through parents with more than 1 metabolite
  {
    zz <- Brnch[[pname]]
    for (brn in names(zz))                      # loop through the metabolites of this parent
    {
      iB <- iB + 1                              # simple counter to move to next element of Bstart/stop one at a time
      Bstart[iB] <- iF                          # same start value for the set of metabolites for each parent
      Bstop[iB] <- iF + length(zz[[brn]]) - 1   # stop value increases for each metabolite of the parent
      iF <- iF + length(zz[[brn]])
    }
  }

  F <- matrix(0.0, nrow=Ndelta, ncol=Ndelta,    # row and columns are all parent-metabolie pairs
              dimnames=list(nms, nms))          # with names corresp to the indices in Xref
  iF <- 1
  for(pname in names(Brnch))                    # loop through parents
  {
    i <- match(pname, rownames(Xref))           # row in Xref corresponding to parent
    zz <- Brnch[[pname]]
    for (brn in names(zz))                      # loop through metabolites of parent
    {
      J <- match(zz[[brn]], colnames(Xref))     # column indices in Xref corresponding to metabolites
      Fcol <- paste(i,J,sep=":")                # parent row : metabolite column for all pairs (e.g. parent 1 metabolites 5,6 -> c(1:5,1:6))
      for (k in 1:length(J))                    # loop over columns in Xref for metabolites
      {
        K <- k + iF - 1                         # K is not used.  Think it is just carried over from copy/past similar function above
        F[Fcol[k], Fcol] <- 1                   #
      }
      iF <- iF + length(J)
    }
  }

  Phi <- Xref

  ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## JAGS inputs: nhanesdata and inits
  ##   -44 rows of NewMeasured have PosteriorShape = "normal", so
  ##   the below fits are just focusing on these chemicals, which
  ##   determines mixmu, mixtau, and mixpi. ly for non-normal seems
  ##   to be estimated (# of samples - # not measured) but don't always
  ##   trend with the loggm column values.  Values SS, lod are non-normal
  ##   only.
  ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (doONEchain)
  {
    nhanesdata <- list()
    ## Now, estimate the mixture distribution for
    lsdlog <- Measured$lsdlog[1:Mn]
    ## Fit 1 and 2 component mixtures, pick 1 unless 2 is significantly better.
    fit1 <- norMixMLE(lsdlog, 1)
    fit2 <- norMixMLE(lsdlog, 2)
    X2 <- 2*(attr(fit2, "loglik") - attr(fit1, "loglik"))
    Pval <- pchisq(X2,3,lower.tail=FALSE)
    Est <- if (Pval < 0.05) fit2 else fit2
    nhanesdata$mixmu <- unname(Est[,"mu"])
    nhanesdata$mixtau <- unname(1/(Est[,"sigma"]^2))
    nhanesdata$mixpi <- unname(Est[,"w"])

    ## ly - the log of the observed metabolite quantiles
    nhanesdata$ly <-
      ly <- unname(c(Measured$loggm[1:Mn],(Measured$Sample_Size[(Mn+1):M] - Measured$BelowLOD[(Mn+1):M])))
    ## se - precision of measured metabolite values, from reported confidence limits on quantiles
    nhanesdata$SS <- unname(Measured$Sample_Size[(Mn+1):M])
    nhanesdata$lod <- log(unname(Measured$LOD[(Mn+1):M]))
    nhanesdata$Nn <- nrow(Est)
    nhanesdata$se <- Measured$loggm_se
    ## Some counts:
    ## N is the number of parents
    nhanesdata$N <- nrow(pred.data)
    ## M is the number of metabolites
    nhanesdata$M <- M
    nhanesdata$Mn <- Mn
    ## indx - indx[i,1:2] is the row and column in Phi for the ith delta
    nhanesdata$indx <- index
    nhanesdata$NBranches <- NBranches
    nhanesdata$Bstart <- Bstart
    nhanesdata$Bstop <- Bstop
    nhanesdata$Alpha <- rep(1, Ndelta)

    ## Phi
    nhanesdata$Phi <- Phi
    nhanesdata$Ndelta <- Ndelta
    ##  nhanesdata$F <- F

    ## ---------------------------------------------------------------------
    ## Initialization
    inits <- list(list())
    ## inits$lP: get estimates for lU, then work backwards to lP
    ## assuming equal distribution of parent to metabolites when unknown
    lU <- nhanesdata$ly
    ## But the last M - Mn values are really counts above lod.  Assuming an lsd at the
    ## middle of the distribution of lsdlog, estimate the value of lU consistent with being above
    ## lod.
    if (Mn < M)
    {
      sd <- exp(median(lsdlog))
      ## This doesn't always work.  Why?
      for (j in (Mn + 1):M)
      {
        lU[j] <- nhanesdata$lod[j - Mn] - qnorm(1 - max(Measured$Prabove[j],.Machine$double.eps)) * sd
      }
    }
    XX <- Phi
    for (i in 1:nrow(XX))
    {
      nNA <- sum(is.na(XX[i,]))
      if (nNA > 0)
      {
        XX[i,is.na(XX[i,])] <- 1/nNA
      }
    }
    XXt <- t(XX)

    obfun <- function(p)
    {
      sum((log(XXt %*% exp(p)) - lU)^2)
    }

    pstart <- as.vector(ginv(XXt) %*% exp(lU))
    pstart[pstart< 0] <- .1
    pstart <- log(pstart)
    out <- optim(pstart, obfun, control=list(maxit=500000))

    inits[[1]]$lP <- out$par

    inits[[1]]$sd.dum <- sd(inits[[1]]$lP)
    phi <- numeric(nhanesdata$Ndelta)
    for (i in 1:NBranches)
    {
      phi[nhanesdata$Bstart[i]:nhanesdata$Bstop[i]] <-
        1/(nhanesdata$Bstop[i] - nhanesdata$Bstart[i] + 1)
    }
    inits[[1]]$phi <- phi
    inits[[1]]$lPmu <- mean(inits[[1]]$lP)

    doinits <- function(chain) {inits[[chain]]}
    assign("inits", inits, envir=environment(doinits))

    ## -------------------------------------------------
    ## Run JAGS

    ## Now, run for 50000 samples, thinned every 100, and dump lP, lU, and delta
    ## as well.  We will use these to get starting values for three new chains.
    if (DEBUGGING)
    {
      N.Iters <- 2000
      N.Thin <- 1
      N.Burnin <- 1000
    } else {
      N.Iters <- 1000
      N.Thin <- 10
      N.Burnin <- 1000
    }

    ## Test code
    phi <- inits[[1]]$phi
    Phi <- nhanesdata$Phi
    indx <- nhanesdata$indx
    for (i in 1:Ndelta) {
      Phi[indx[i,1],indx[i,2]] <- phi[i]
    }

    model.jags <- jags.model(bayes_model, data=nhanesdata, inits=doinits,
                             n.chains=1,n.adapt=N.Burnin)

    out.samples <- coda.samples(model.jags,
                                variable.names=c("lP","phi","lU","lPmu","sd.V"),
                                n.iter=N.Iters * N.Thin, thin=N.Thin)

    nms <- varnames(out.samples)
    #save(out.samples, nhanesdata, inits, file=paste("OnlyPparms1_", SUBPOP, ".RData", sep = ""))
  }
  ## --------------------------------------------------------------------------------
  ## This is the initial, hopefully converged, chain.  Now, use these to start with
  ## and run three new chains.

  print("Set up initial conditions")

  ## getInits finds extreme values
  getInits <- function(X, ninits=3, P=0.05)
  {
    ## X is an mcmc object, essentially a matrix.
    dta <- scale(X)
    ## Ultimately, we are looking for ninits indices into the elements of X
    outinits <- numeric(ninits)
    ## Pick ninits random directions.  Project the data in dta on each direction
    ## and pick the Pth quantile (since the original direction u was random
    ## u and -u are equally likely).
    for (i in 1:ninits)
    {
      u <- rnorm(ncol(dta))
      u <- u/sqrt(sum(u*u))
      lngt <- dta %*% u
      indx <- seq(along=lngt)
      keep <- round(P*length(indx)) ## this works fine as long as length(indx) > 20
      oind <- order(lngt)
      outinits[i] <- oind[keep]
    }
    outinits
  }
  ## Get 3 hyperdispersed points.
  out.coda <- out.samples[[1]]

  Indx <- getInits(out.coda,3)

  getsamps <- function(nm)
  {
    nm <- gsub("\\.","\\\\.",nm)
    cindx <- grep(nm, colnames(out.coda))
    out.coda[,cindx,drop=FALSE]
  }

  inits2 <- list(list(),list(),list())
  for (i in 1:3)
  {
    inits2[[i]]$lP <- getsamps("^lP\\[")[Indx[i],]
    inits2[[i]]$sd.dum <- getsamps("^sd.V")[Indx[i],]
    inits2[[i]]$phi <- getsamps("phi")[Indx[i],]
    inits2[[i]]$lPmu <- getsamps("lPmu")[Indx[i],]
    inits2[[i]]$.RNG.name="lecuyer::RngStream"
    inits2[[i]]$.RNG.seed=randomNumbers(n=1, min=1, max=1e+06, col=1)
  }

  doinits <- function(chain){inits2[[chain]]}
  env1 <- new.env()
  assign("inits2", inits2, pos=env1)
  environment(doinits) <- env1

  if (DEBUGGING)
  {
    NIters <- 2000
    Thin <- 1
    NBurnin <- 1000
  } else {
    NIters <- 1000
    Thin <- 500
    NBurnin <- 5000
  }

  print("Start main computation")
  out.samps3R <- foreach(i = 1:3) %dopar% {
    model <- jags.model(bayes_model, data=nhanesdata,
                        inits=inits2[i], n.chains=1, n.adapt=NBurnin)
    result <- coda.samples(model, variable.names=c("lP","phi","lU","lPmu","sd.V"),
                           n.iter=NIters*Thin, thin=Thin)
    return(result)
  }

  out.coda3R <- do.call("mcmc.list", lapply(out.samps3R, function(z) z[[1]]))

  print(paste("Save final result in the file OnlyPparms3_", SUBPOP, "_",
              format(Sys.time(), "%Y-%m-%d"), ".RData", sep = ""))
  save(out.samps3R, out.coda3R, file=paste("OnlyPparms3_", SUBPOP, "_",
                                           format(Sys.time(), "%Y-%m-%d"), ".RData", sep = ""))
  return(out.samps3R, out.coda3R)

}


