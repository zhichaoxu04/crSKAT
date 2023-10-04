#' genData.R
#'
#' Generate competing risks interval-censored data under the proportional odds given a baseline hazard function and
#' some information about observation times.
#'
#' @param seed Random seed.
#' @param n The number of observations to generate (sample size).
#' @param alpha1 1st parameter for Gompertz model.
#' @param alpha2 2nd parameter for Gompertz model.
#' @param beta1 3rd parameter for Gompertz model.
#' @param beta2 4th parameter for Gompertz model.
#' @param obsTimes Vector of the intended observation times.
#' @param probMiss The probability of missing any given visit (default is 0.1).
#' @param windowHalf The amount of time before or after the intended obsTimes that a visit might take place (default is 1).
#'
#' @return A list with the elements:
#' \item{tempTime}{A vector of exact event times.}
#' \item{tempType}{A vector of exact event types.}
#' \item{deltaVecSimple}{A vector of whether the event was observed after follow-up started (t>0).}
#' \item{deltaVec}{A vector of the event cause (cause 1, cause 2, or right-censored 0).}
#' \item{leftTimes}{A vector of left side of interval times.}
#' \item{rightTimes}{A vector of right side of interval times.}
#'
#' @importFrom stats runif
#' @importFrom stats rbinom
#'
#' @export
#' @examples
#' seed <- 2022
#' n <- 1000
#' alpha1 <- -0.058
#' alpha2 <- -0.035
#' beta1 <- 0.03
#' beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2
#' outcomeDat <- genData(seed=seed, n=n, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2,
#' obsTimes=seq(from=4, to=28, by=4), probMiss=0.1, windowHalf=1)
#'
genData <- function(seed=NULL, n, alpha1, alpha2, beta1, beta2,
                    obsTimes=seq(from=4, to=28, by=4), probMiss=0.1, windowHalf=1) {
  if (!is.null(seed)) {
    set.seed(seed)
  } else{
    set.seed(sample(1:10000, 1))
  }

  # assume as in Hudgens et al. paper that \sum Fk(\infty) = 1, so one of events will always happen,
  # of course alphas and betas are the determinants of whether this is true
  pi1 <- 1 - exp(beta1 / alpha1)
  pi2 <- 1 - exp(beta2 / alpha2)
  tempType <- stats::rbinom(n=n, size=1, prob=pi2) + 1
  tempUnif <- stats::runif(n=n)
  tempTime <- ifelse(tempType == 1, log( 1 - log((1 - tempUnif * pi1)) * alpha1 / beta1 ) / alpha1,
                     log( 1 - log((1 - tempUnif * pi2)) * alpha2 / beta2 ) / alpha2)

  # observation time
  # simulate missing visit
  nVisits <- length(obsTimes)
  madeVisit <- matrix(data=stats::rbinom(n=n*nVisits, size=1, prob=(1 - probMiss)), nrow=n, ncol=nVisits)

  # make sure there is at least one visit for each subject
  nMadeVisits <- apply(madeVisit, 1, sum)
  zeroVisits <- which(nMadeVisits == 0)
  while (length(zeroVisits) > 0) {
    madeVisit[zeroVisits, ] <- matrix(data = stats::rbinom(n=length(zeroVisits) * nVisits, size=1,
                                                           prob=(1 - probMiss)), nrow=length(zeroVisits), ncol=nVisits)
    nMadeVisits <- apply(madeVisit, 1, sum)
    zeroVisits <- which(nMadeVisits == 0)
  }

  visitTime <- sweep(matrix(data=stats::runif(n=n*nVisits, min=-windowHalf, max=windowHalf), nrow=n, ncol=nVisits),
                     MARGIN=2, STATS=obsTimes, FUN="+")
  # get all visits for each subject
  allVisits <- madeVisit * visitTime
  # make the interval for each subject - USING creatIntNew!
  allInts <- t(mapply(FUN=createIntNew, obsTimes = data.frame(t(allVisits)), eventTime=tempTime))
  leftTimes <- allInts[, 1]
  rightTimes <- allInts[, 2]
  # new indicator of whether event was observed
  deltaVecSimple <- ifelse(rightTimes > tempTime, 1, 0)
  deltaVec <- deltaVecSimple * tempType

  return(list(tempTime = tempTime, tempType = tempType, deltaVecSimple = deltaVecSimple, deltaVec = deltaVec,
              leftTimes = leftTimes, rightTimes = rightTimes))
}
