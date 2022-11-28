#' createIntNew.R
#'
#' Turn the actual observed outcome times and observation times into interval-censored
#' outcomes for each subject. Apply this with mapply over a data.frame of visit times, pass in the exact times.
#'
#' @param obsTimes A vector of all the times a subject is observed.
#' @param eventTime The exact event time for the subject.
#'
#' @return A 2*1 vector which is the interval of the event time.
#'
#'
#' @export
#' @examples
#' obsTimes <- 1:10
#' eventTime <- 4.4
#' createIntNew(obsTimes, eventTime)
#'
createIntNew <- function(obsTimes, eventTime) {
  # order the times in case the random portion causes them to go out of order
  orderedTimes <- sort(obsTimes)
  # left end of interval
  minIdx <- which(orderedTimes < eventTime)
  if (length(minIdx) == 0) {
    minTime <- 0
  } else {
    minTime <- orderedTimes[max(minIdx)]
  }
  # right end of interval
  maxIdx <- which(orderedTimes >= eventTime)
  if (length(maxIdx) == 0) {
    maxTime <- max(orderedTimes)  # Right Censored
    minTime <- 0
  } else {
    maxTime <- orderedTimes[min(maxIdx)]
  }

  return(c(minTime, maxTime))
}
