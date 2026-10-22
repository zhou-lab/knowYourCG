#' Master data frame for all object to cache
#'
#' This is an internal object, updated on every release from the
#' knowYourCG records on ExperimentHub.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{name}{Character. Resource name/title.}
#'   \item{eh}{Character. ExperimentHub record ID.}
#' }
#' 
#' @name df_master
#' @docType data
#' @keywords internal
#' @return master sheet of knowYourCG objects
NULL

cacheEnv <- new.env()
