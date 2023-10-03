## code to prepare `exampleData` dataset goes here

genDataraw <- function(seed=NULL, n, alpha1, alpha2, beta1, beta2) {
  if (!is.null(seed)) {
    set.seed(seed)
  }else {
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

  return(list(tempTime = tempTime, tempType = tempType))
}

alpha1 <- -0.058
alpha2 <- -0.035
beta1 <- 0.03
beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2

exampleData <- genDataraw(seed=NULL, n=1000, alpha1, alpha2, beta1, beta2)
usethis::use_data(exampleData, overwrite = TRUE)
