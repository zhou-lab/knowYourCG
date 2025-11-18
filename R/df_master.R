#' Master data frame for all object to cache
#'
#' This is an internal object which will be updated on every new release
#' library(ExperimentHub)
#' eh <- query(ExperimentHub(localHub=FALSE), "knowYourCG")
#' eh <- query(ExperimentHub(localHub=FALSE), "sesameData") # older data
#' data.frame(name=eh$title, eh=names(eh))
#'
#' Cache location is default to
#' /Users/zhouw3/Library/Caches/org.R-project.R/R/ExperimentHub/
#'
#' @name df_master
#' @docType data
#' @return master sheet of knowYourCG objects
NULL

cacheEnv <- new.env()
