#' scoreEqSpline.R
#'
#' Calculate score function using cubic spline for baseline cumulative hazard.
#'
#' @param x (p+nknots+2)*1 vector of fitted coefficients under null model.
#' @param leftDmat n*(p+nknots+2) design matrix for left end of interval.
#' @param rightDmat n*(p+nknots+2) design matrix for right end of interval.
#' @param gSummed n*1 vector of summed genotype matrix (default is NULL).
#' @param gMat n*q genotype matrix.
#' @param estG Whether to include genotype matrix into design matrix (default is FALSE).
#' @param leftTimes n*1 vector of left side of interval times.
#' @param deltaVec n*1 vector of the event cause (cause 1 or 2; or right-censored 0).
#'
#' @return A score vector
#'
#'
#' @export
#' @examples
#' seed <- 2022
#' n <- 1000
#' q <- 20
#' alpha1 <- -0.058
#' alpha2 <- -0.035
#' beta1 <- 0.03
#' beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2
#' outcomeDat <- genData(seed=seed, n=n, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)
#' xMat <- matrix(data=rnorm(n), nrow=n)
#' gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.3), ncol=q)
#' gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
#' obsTimes <- 1:5
#' lt <- outcomeDat$leftTimes
#' rt <- outcomeDat$rightTimes
#' obsInd <- outcomeDat$deltaVecSimple
#' dmats <- makeICdmat(xMat = xMat, lt = lt, rt = rt, obsInd = obsInd, quant_r = NULL, nKnots = 1)
#' forUgamma <- scoreEqSpline(x=runif(9, 0.01, 0.02), leftDmat = dmats$left_dmat,
#' rightDmat = dmats$right_dmat, leftTimes = lt, deltaVec = outcomeDat$deltaVec,
#' gSummed = gSummed, gMat = gMat, estG = FALSE)
#'
scoreEqSpline <- function(x, leftDmat, rightDmat, gSummed=NULL, gMat=NULL, estG=FALSE, leftTimes, deltaVec) {

  # different designs for different outcomes
  # here there should be no "infinity" type terms in the rightDmat

  leftDmat2 <- cbind(leftDmat, gSummed)
  rightDmat2 <- cbind(rightDmat, gSummed)
  if (!is.null(gMat) & estG) {
    leftDmat1 <- cbind(leftDmat, gMat)
    rightDmat1 <- cbind(rightDmat, gMat)
  } else {
    leftDmat1 <- leftDmat
    rightDmat1 <- rightDmat
  }

  theta1 <- x[1:ncol(leftDmat1)]
  theta2 <- x[(ncol(leftDmat1)+1):length(x)]

  # separate deltas
  delta0 <- ifelse(deltaVec == 0, 1, 0)
  delta1 <- ifelse(deltaVec == 1, 1, 0)
  delta2 <- ifelse(deltaVec == 2, 1, 0)

  # H and S terms
  H1L <- ifelse(leftTimes == 0, 0, exp(leftDmat1 %*% theta1))
  #H1L <- exp(leftDmat %*% theta1)
  # have to use 999 because Inf * 0  = Inf
  #H1R <- ifelse(deltaVec == 0, 999, exp(rightDmat %*% theta1))
  H1R <- exp(rightDmat1 %*% theta1)
  S1L <- ifelse(leftTimes == 0, 1, exp(-H1L))
  S1R <- exp(-H1R)
  H2L <- ifelse(leftTimes == 0, 0, exp(leftDmat2 %*% theta2))
  #H2L <- exp(leftDmat %*% theta2)
  #H2R <- ifelse(deltaVec == 0, 999, exp(rightDmat %*% theta2))
  H2R <- exp(rightDmat2 %*% theta2)
  S2L <- ifelse(leftTimes == 0, 1, exp(-H2L))
  S2R <- exp(-H2R)

  # to be swept
  L1term <- ifelse(leftTimes == 0, 0, -delta1 * H1L * S1L / (S1L - S1R))
  R1term <- delta1 * H1R * S1R / (S1L - S1R)
  E1term <- -delta0 * H1R * S1R / (S1R + S2R - 1)
  L2term <- ifelse(leftTimes == 0, 0, -delta2 * H2L * S2L / (S2L - S2R))
  R2term <- delta2 * H2R * S2R / (S2L - S2R)
  E2term <- -delta0 * H2R * S2R / (S1R + S2R - 1)
  # sweep it
  U1swept <- sweep(leftDmat1, MARGIN=1, STATS=L1term, FUN="*") +
    sweep(rightDmat1, MARGIN=1, STATS=R1term, FUN="*") +
    sweep(rightDmat1, MARGIN=1, STATS=E1term, FUN="*")
  U1 <- colSums(U1swept)
  U2swept <- sweep(leftDmat2, MARGIN=1, STATS=L2term, FUN="*") +
    sweep(rightDmat2, MARGIN=1, STATS=R2term, FUN="*") +
    sweep(rightDmat2, MARGIN=1, STATS=E2term, FUN="*")
  U2 <- colSums(U2swept)

  # this is just to get the score equations for \gamma under the null
  if (!is.null(gMat) & !estG) {
    UgSwept <- sweep(gMat, MARGIN=1, STATS=L1term + R1term + E1term, FUN="*")
    uVec <- c(U1, colSums(UgSwept), U2)
  } else {
    uVec <- c(U1, U2)
  }
  return(uVec)

}
