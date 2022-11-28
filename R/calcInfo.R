#' calcInfo.R
#'
#' Calculate Fisher information matrix.
#'
#' @param leftDmat n*(p+nknots+2) design matrix for left end of interval.
#' @param rightDmat n*(p+nknots+2) design matrix for right end of interval.
#' @param gSummed n*1 vector of summed genotype matrix.
#' @param gMat n*q genotype matrix (default is NULL).
#' @param leftTimes n*1 vector of left side of interval times.
#' @param theta1 n*1 vector of fitted coefficient for cause 1.
#' @param theta2 n*1 vector of fitted coefficient for cause 2.
#' @param deltaVec n*1 vector of the event cause (cause 1 or 2; or right-censored 0).
#'
#' @return A Fisher information matrix for model coefficients
#'
#'
#' @export
#' @examples
#' seed <- 2022
#' n <- 3000
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
#' iMatPartial <- calcInfo(leftDmat = dmats$left_dmat, rightDmat = dmats$right_dmat, leftTimes = lt,
#' theta1 = runif(4, 0.01, 0.02), theta2 = runif(5, 0.01, 0.02), deltaVec = outcomeDat$deltaVec,
#' gSummed=gSummed, gMat = NULL)
#'
calcInfo <- function(leftDmat, rightDmat, gSummed, gMat = NULL, leftTimes, theta1, theta2, deltaVec) {

  # design matrices
  leftDmat2 <- cbind(leftDmat, gSummed)
  rightDmat2 <- cbind(rightDmat, gSummed)
  if (!is.null(gMat)) {
    leftDmat1 <- cbind(leftDmat, gMat)
    rightDmat1 <- cbind(rightDmat, gMat)
  } else {
    leftDmat1 <- leftDmat
    rightDmat1 <- rightDmat
  }

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

  # dt1dt1 part
  d11ToSweep1 <- ifelse(leftTimes == 0, 0, delta1 * H1L * S1L * (1 - H1L)) / (S1L - S1R)
  d11ToSweep2 <- delta1 * H1R * S1R * (1 - H1R) / (S1L - S1R)
  d11ToSweep3 <- ifelse(leftTimes == 0, 0, delta1 * H1L * S1L) / (S1L - S1R)^1
  d11ToSweep4 <- delta1 * H1R * S1R / (S1L - S1R)^1
  d11ToSweep5 <- delta0 * H1R * S1R * (1 - H1R) / (S1R + S2R - 1)
  d11ToSweep6 <- delta0 * H1R * S1R / (S1R + S2R - 1)^1

  dt1dt1 <- -t(leftDmat1) %*% sweep(leftDmat1, MARGIN=1, STATS=d11ToSweep1, FUN="*") +
    t(rightDmat1) %*% sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep2, FUN="*") +
    -t(sweep(-leftDmat1, MARGIN=1, STATS=d11ToSweep3, FUN="*") + sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep4, FUN="*")) %*%
    (sweep(-leftDmat1, MARGIN=1, STATS=d11ToSweep3, FUN="*") + sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep4, FUN="*")) +
    -t(rightDmat1) %*% sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep5, FUN="*") +
    -t(sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep6, FUN="*")) %*% sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep6, FUN="*")

  # dt2dt2 part
  d22ToSweep1 <- ifelse(leftTimes == 0, 0, delta2 * H2L * S2L * (1 - H2L)) / (S2L - S2R)
  d22ToSweep2 <- delta2 * H2R * S2R * (1 - H2R) / (S2L - S2R)
  d22ToSweep3 <- ifelse(leftTimes == 0, 0, delta2 * H2L * S2L) / (S2L - S2R)^1
  d22ToSweep4 <- delta2 * H2R * S2R / (S2L - S2R)^1
  d22ToSweep5 <- delta0 * H2R * S2R * (1 - H2R) / (S1R + S2R - 1)
  d22ToSweep6 <- delta0 * H2R * S2R / (S1R + S2R - 1)^1
  dt2dt2 <- -t(leftDmat2) %*% sweep(leftDmat2, MARGIN=1, STATS=d22ToSweep1, FUN="*") +
    t(rightDmat2) %*% sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep2, FUN="*") +
    -t(sweep(-leftDmat2, MARGIN=1, STATS=d22ToSweep3, FUN="*") + sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep4, FUN="*")) %*%
    (sweep(-leftDmat2, MARGIN=1, STATS=d22ToSweep3, FUN="*") + sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep4, FUN="*")) +
    -t(rightDmat2) %*% sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep5, FUN="*") +
    -t(sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep6, FUN="*")) %*% sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep6, FUN="*")

  # off diagonal part
  dt1dt2 <- -t(sweep(rightDmat1, MARGIN=1, STATS=delta0 * H1R * S1R / (S1R + S2R - 1)^1, FUN="*")) %*%
    sweep(rightDmat2, MARGIN=1, STATS=delta0 * H2R * S2R / (S1R + S2R - 1)^1, FUN="*")

  # put it together
  iMat <- rbind(cbind(dt1dt1, dt1dt2), cbind(t(dt1dt2), dt2dt2))

  return(iMat)
}
