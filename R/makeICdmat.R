#' makeICdmat.R
#'
#' Puts together the entire design matrix for both the left and right ends of the
#' interval, pasting together the non-genetic covariates with the cubic spline basis.
#'
#' @param xMat n*p matrix of non-genetic covariates.
#' @param lt n*1 vector with left end of intervals (min is 0).
#' @param rt n*1 vector with right end of intervals.
#' @param obsInd n*1 vector of whether the event was observed before last follow-up.
#' @param quant_r Quantiles of time to use in constructing the spline, pass in if doing bootstrap.
#' @param nKnots Number of knots to use for cubic spline basis (default is 1).
#'
#' @return A list with the elements:
#' \item{right_dmat}{n*(p+nKnots+2) design matrix for right end of interval.}
#' \item{left_dmat}{n*(p+nKnots+2) design matrix for left end of interval.}
#' \item{quant_r}{Quantiles used for constructing spline.}
#'
#' @importFrom stats quantile
#'
#' @export
#' @examples
#' seed <- 2022
#' n <- 1000
#' alpha1 <- -0.058
#' alpha2 <- -0.035
#' beta1 <- 0.03
#' beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2
#' outcomeDat <- genData(seed=seed, n=n, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
#' xMat <- matrix(data=rnorm(n), nrow=n)
#' obsTimes <- 1:5
#' lt <- outcomeDat$leftTimes
#' rt <- outcomeDat$rightTimes
#' obsInd <- outcomeDat$deltaVecSimple
#' makeICdmat(xMat = xMat, lt = lt, rt = rt, obsInd = obsInd, quant_r = NULL, nKnots = 1)
#'
makeICdmat <- function(xMat, lt, rt, obsInd, quant_r = NULL, nKnots = 1) {

  # place the knots at equally spaced quantiles
  if (is.null(quant_r)) {
    quant_r <- stats::quantile(log(rt[obsInd == 1]), probs=seq(from=0, to=1, length.out=nKnots+2))
  }

  # a0 and a1 are always there
  right_a0 <- 1
  #right_a1 <- ifelse(obsInd == 0, 999, log(rt))
  right_a1 <- log(rt)
  if (is.null(xMat)) {
    right_dmat <- cbind(right_a0, right_a1)
  } else {
    right_dmat <- cbind(xMat, right_a0, right_a1)
  }

  # if lt = 0, then the cumulative hazard is necessarily 0 so set it all to 0
  left_a0 <- ifelse(lt == 0, 0, 1)
  left_a1 <- ifelse(lt == 0, 0, log(lt))
  if (is.null(xMat)) {
    left_dmat <- cbind(left_a0, left_a1)
  } else {
    left_dmat <- cbind(xMat, left_a0, left_a1)
  }

  kmax <- max(quant_r)
  kmin <- min(quant_r)
  # place the knots
  for (j in 1:nKnots) {
    ej <- (kmax - quant_r[j+1]) / (kmax - kmin)

    #right_aj <-  ifelse(obsInd == 0, 999,
    #                    pmax(0, (right_a1 - quant_r[j+1])**3) - ej * pmax(0, (right_a1 - kmin)**3) -
    #                      (1 - ej) * pmax(0, (right_a1 - kmax)**3))
    right_aj <- pmax(0, (right_a1 - quant_r[j+1])**3) - ej * pmax(0, (right_a1 - kmin)**3) -
      (1 - ej) * pmax(0, (right_a1 - kmax)**3)
    right_dmat <- cbind(right_dmat, right_aj)

    left_aj <- ifelse(lt == 0, 0, pmax(0, (left_a1 - quant_r[j+1])**3) - ej * pmax(0, (left_a1 - kmin)**3) -
                        (1 - ej) * pmax(0, (left_a1 - kmax)**3))
    left_dmat <- cbind(left_dmat, left_aj)
  }

  return(list(right_dmat=right_dmat, left_dmat=left_dmat, quant_r = quant_r))
}
