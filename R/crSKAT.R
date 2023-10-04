#' crSKAT.R
#'
#' Calculate the test statistic and p-value for interval-censored competing risks SKAT.
#'
#' @param leftDmat n*(p+nknots+2) design matrix for left end of interval.
#' @param rightDmat n*(p+nknots+2) design matrix for right end of interval.
#' @param leftTimes n*1 vector of left side of interval times.
#' @param deltaVec n*1 vector of the event cause (cause 1 or 2; or right-censored 0).
#' @param gMat n*q genotype matrix.
#' @param gSummed n*1 vector of summed genotype matrix.
#' @param null_beta (p+nknots+2)*1 vector of coefficients for null model.
#' @param pvalue If TRUE, then find the p-value
#'
#' @return A list with the elements:
#' \item{p_SKAT}{ICSKAT p-value}
#' \item{p_burden}{IC burden test p-value}
#' \item{complex}{Indicator of whether the SKAT variance matrix was positive definite}
#' \item{sig_mat}{The covariance matrix of the score equations for genetic effects when treated as fixed effects}
#' \item{skatQ}{SKAT test statistic}
#' \item{burdenQ}{Burden test statistic}
#' \item{Ugamma}{Score vector}
#' \item{lambdaQ}{Vector of eigenvalues of variance matrix}
#' \item{null_beta}{The fitted null parameters}
#' \item{err}{Will be 0 for no error, 22 if had to adjust parameters on CompQuadForm (totally normal), or 99 if NA in variance matrix. ICSKATwrapper will return 1 here if the null fit has an error}
#' \item{errMsg}{Explains error code, blank string if no error}
#'
#' @importFrom CompQuadForm davies
#' @importFrom stats pchisq
#'
#' @export
#' @examples
#' n <- 30000
#' q <- 500
#' alpha1 <- -0.1
#' alpha2 <- -0.3
#' beta1 <- 0.03
#' beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2
#' outcomeDat <- genData(seed=NULL, n=n, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
#' xMat <- matrix(data=rnorm(n), nrow=n)
#' gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.3), ncol=q)
#' gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
#' obsTimes <- 3:7
#' lt <- outcomeDat$leftTimes
#' rt <- outcomeDat$rightTimes
#' obsInd <- outcomeDat$deltaVecSimple
#' deltaVec <- outcomeDat$deltaVec
#' dmats <- makeICdmat(xMat = xMat, lt = lt, rt = rt, obsInd = obsInd, quant_r = NULL, nKnots = 1)
#' leftDmat <- dmats$left_dmat
#' rightDmat <- dmats$right_dmat
#' nullFit <- crSKAT_fit_null(init_beta=rep(0, 9), leftDmat=leftDmat, rightDmat=rightDmat,
#' deltaVec=deltaVec, leftTimes=lt, gSummed=gSummed, method="Broyden")
#' out <- crSKAT(leftDmat=dmats$left_dmat, rightDmat=rightDmat, leftTimes=lt,
#' deltaVec=deltaVec, gMat=gMat, gSummed=gSummed, null_beta=nullFit$beta_fit, pvalue=TRUE)
#'
crSKAT <- function(leftDmat, rightDmat, leftTimes, deltaVec, gMat, gSummed, null_beta, pvalue=TRUE) {
  q <- ncol(gMat)
  numCov <- ncol(leftDmat)

  if(is.na(null_beta)[1]){
    return(list(lambdaQ=NA, p_SKAT=NA, p_burden=NA, skatQ=NA, Ugamma=NA, null_beta = null_beta,
                burdenQ=NA, sig_mat = NA, complex=NA, err=3, errMsg="Null beta input is NA"))
  }

  forUgamma <- scoreEqSpline(x=null_beta, leftDmat=leftDmat,rightDmat=rightDmat, leftTimes=leftTimes,
                             deltaVec=deltaVec, gSummed=gSummed, gMat=gMat, estG=FALSE)

  iMatPartial <- calcInfo(leftDmat=leftDmat, rightDmat=rightDmat, leftTimes=leftTimes,
                          theta1=null_beta[1:numCov],
                          theta2=null_beta[(numCov+1):length(null_beta)],
                          deltaVec=deltaVec, gSummed=gSummed, gMat=NULL)


  forIgt <- calcInfo(leftDmat = leftDmat, rightDmat = rightDmat, leftTimes = leftTimes,
                     theta1 = c(null_beta[1:numCov], rep(0, q)),
                     theta2 = null_beta[(numCov+1):length(null_beta)],
                     deltaVec = deltaVec, gSummed = gSummed, gMat = gMat)

  skatQ <- t(forUgamma[(numCov+1):(numCov+q)]) %*% forUgamma[(numCov+1):(numCov+q)]
  burdenQ <- (sum(forUgamma[(numCov+1):(numCov+q)]))^2
  Itt <- -iMatPartial
  Igg <- -forIgt[(numCov+1):(numCov+q), (numCov+1):(numCov+q)]
  Igt <- -forIgt[c(1:numCov, ((numCov+1)+q):nrow(forIgt)), (numCov+1):(numCov+q)]
  if (length(which(is.na(Igg))) > 0 | length(which(is.na(Igt))) > 0 | length(which(is.na(Itt))) > 0) {
    return(list(lambdaQ=NA, p_SKAT=NA, p_burden=NA, skatQ=NA, Ugamma=forUgamma, null_beta = null_beta,
                burdenQ=NA, sig_mat = NA, complex=NA, err=99, errMsg="NA in variance matrix"))
  }

  sig_mat <- Igg - t(Igt) %*% solve(Itt) %*% (Igt)

  if (length(which(is.na(sig_mat))) > 0) {
    return(list(lambdaQ=NA, p_SKAT=NA, p_burden=NA, skatQ=NA, Ugamma=forUgamma, null_beta = null_beta,
                burdenQ=NA, sig_mat = NA, complex=NA, err=99, errMsg="NA in variance matrix"))
  }

  errCode <- 0
  errMsg <- ""

  # calculate p-value
  if (pvalue) {
    lambdaQ <- eigen(sig_mat)$values
    p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
    # as noted in the CompQuadForm documentation, sometimes you need to play with acc or lim parameters
    # to get a p-value between 0 and 1
    if (!is.na(p_SKAT)) {
      if (p_SKAT > 1) {
        paramDF <- data.frame(expand.grid(lim = c(10000, 20000, 50000), acc=c(1e-7, 1e-6, 1e-5, 1e-4)))
        paramCounter <- 1
        while(p_SKAT > 1) {
          tempLim <- paramDF$lim[paramCounter]
          tempAcc <- paramDF$acc[paramCounter]
          p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=tempAcc, lim=tempLim)$Qq
          paramCounter <- paramCounter + 1
          if (paramCounter > nrow(paramDF)) {break}
        }
        errCode <- 22
        errMsg <- "Had to adjust parameters on CompQuadForm"
      }
    }
    B_burden= burdenQ / sum(sig_mat);
    p_burden= 1 - stats::pchisq(B_burden, df = 1)
  } else {
    pSKAT <- NA
    lambdaQ <- 1
    p_burden <- NA
  }

  return(list(lambdaQ=lambdaQ, p_SKAT=p_SKAT, p_burden=p_burden, skatQ=skatQ, Ugamma=forUgamma, null_beta = null_beta,
              burdenQ=burdenQ, sig_mat = sig_mat, complex=is.complex(lambdaQ), err=errCode, errMsg=errMsg))

}
