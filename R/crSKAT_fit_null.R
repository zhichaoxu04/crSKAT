#' crSKAT_fit_null.R
#'
#' Fit the null model (cubic basis spline for baseline cumulative hazard and coefficients
#' for non-genetic coefficiens) for interval-censored competing risks SKAT.
#'
#' @param init_beta (p+nknots+2)*1 vector of coefficients to initialize.
#' @param leftDmat n*(p+nknots+2) design matrix for left end of interval.
#' @param rightDmat n*(p+nknots+2) design matrix for right end of interval.
#' @param deltaVec n*1 vector of the event cause (cause 1 or 2; or right-censored 0).
#' @param leftTimes n*1 vector of left interval times.
#' @param gSummed n*1 vector of summed genotype matrix.
#' @param allowSingular A logical value indicating if a small correction to the Jacobian when it is singular or too ill-conditioned is allowed.
#' @param method The method to use for finding a solution.
#'
#' @return A list with the elements:
#' \item{beta_fit}{(p+nknots+2)*1 vector of fitted coefficients under null model.}
#' \item{iter}{Number of iterations needed to converge.}
#' \item{times}{Number of times to add variation to the initial coefficients.}
#' \item{Msg}{A string describing the termination code of the nonlinear equations with Broyden or Newton.}
#' \item{termcd}{Termination code as integer.}
#' \item{err}{Value is 0 if no errors and 1 if it can't perform fit.}
#' \item{errMsg}{Empty string if err=0, explains error if there is one.}
#'
#' @importFrom nleqslv nleqslv
#' @importFrom stats rnorm
#'
#' @export
#' @examples
#' seed <- 2022
#' n <- 300000
#' q <- 20
#' alpha1 <- -0.058
#' alpha2 <- -0.035
#' beta1 <- 0.03
#' beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2
#' outcomeDat <- genData(seed=seed, n=n, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
#' deltaVec <- outcomeDat$deltaVec
#' xMat <- matrix(data=rnorm(n), nrow=n)
#' gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.3), ncol=q)
#' gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
#' obsTimes <- 1:5
#' lt <- outcomeDat$leftTimes
#' rt <- outcomeDat$rightTimes
#' obsInd <- outcomeDat$deltaVecSimple
#' dmats <- makeICdmat(xMat = xMat, lt = lt, rt = rt, obsInd = obsInd, quant_r = NULL, nKnots = 1)
#' rightDmat <- dmats$right_dmat
#' crSKAT_fit_null(init_beta=rep(0, 9), leftDmat=dmats$left_dmat, rightDmat=rightDmat,
#' deltaVec=deltaVec, leftTimes=lt, gSummed=gSummed, method="Broyden")
#'
crSKAT_fit_null <- function(init_beta, leftDmat, rightDmat, deltaVec, leftTimes, gSummed, allowSingular=TRUE, method = c("Broyden", "Newton")) {
  if(!method %in% c("Broyden", "Newton")){
    return(list(beta_fit=NA, iter=NA, times=iter, Msg=NA,
                termcd=NA, err=1, errMsg="Specify the method within the list"))
  }

  # solve score equations under null with no random effect
  iter <- 0
  solPartial <- NA
  while(iter < 100 & is.na(solPartial)[1]){
    solPartial <- tryCatch(nleqslv::nleqslv(x = init_beta, fn=scoreEqSpline, leftDmat=leftDmat, rightDmat=rightDmat,
                                            leftTimes=leftTimes, deltaVec=deltaVec, gSummed=gSummed, gMat=NULL, estG=FALSE,
                                            method=method,
                                            control=list(allowSingular=allowSingular)),
                           error = function(e){return(as.numeric(NA))})
    iter <- iter + 1

    # Add some variation to initial value
    set.seed(iter)
    init_beta <- init_beta + rnorm(length(init_beta), 0, 0.01)
  }

  if(is.na(solPartial)[1]){
    return(list(beta_fit=NA, iter=NA, times=iter, Msg=NA,
                termcd=NA, err=1, errMsg="Try different initial values"))
  }

  return(list(beta_fit=as.numeric(solPartial$x), iter=solPartial$iter, times=iter, Msg=solPartial$message,
              termcd=solPartial$termcd, err=0, errMsg=""))

}
