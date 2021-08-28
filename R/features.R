#' calcDatabaseSetStatistics1 calculates features of x
#'
#' @param x Vector of numeric values
#'
#' @return Vector with ~20 different engineered features
#'
#' @import stats
calcDatabaseSetStatistics1 = function(x) {
    a = data.frame(mean=apply(x, 2, mean, na.rm=TRUE),
      median=apply(x, 2, median, na.rm=TRUE),
      var=apply(x, 2, var, na.rm=TRUE),
      sd=apply(x, 2, sd, na.rm=TRUE),
      skew=apply(x, 2, var, na.rm=TRUE),
      iqr=apply(x, 2, IQR, na.rm=TRUE),
      range=apply(x, 2, max, na.rm=TRUE) - apply(x, 2, min, na.rm=TRUE),
      min=apply(x, 2, min, na.rm=TRUE),
      max=apply(x, 2, max, na.rm=TRUE))
  b = apply(x, 2, quantile, na.rm=TRUE, probs=seq(0, 1, 0.1))
  return(cbind(a, t(b)))
}
    

#' calcDatabaseSetStatisticsAll builds dataset for a given filename composed
#' of engineered features from the given databaseSets
#'
#' @param filename String representing a filename to load using readRDS/tbk_data
#' @param databaseSets List of vectors corresponding to probe locations for
#' which the features will be extracted
#'
#' @return Vector for a given sample columns are features across different
#' databaseSets
calcDatabaseSetStatisticsAll = function(betas, databaseSets) {
    a = do.call(cbind, 
            lapply(names(databaseSets),
                   function(databaseSetName) {
                      databaseSet = databaseSets[[databaseSetName]]
                      if (length(databaseSet) >= nrow(betas)) return(FALSE)
                      if (is.numeric(databaseSet)) {
                          probes = names(databaseSet)
                      } else {
                          probes = databaseSet
                      }
                      
                      statistics = suppressWarnings(calcDatabaseSetStatistics1(betas[na.omit(match(probes, rownames(betas))), ]))
                      names(statistics) = unlist(lapply(names(statistics), function(colname) {
                          paste(databaseSetName, colname, sep="-")
                      }))
                      return(statistics)
                  }))
    b = a[, !grepl("FALSE", colnames(a))]
    c = b[, !apply(b, 2, function(x) {any(is.na(x) | is.infinite(x))})]
    return(c)
}


#' buildStatisticDataSet builds dataset for each filename in filenames composed
#' of engineered features from the given databaseSets
#'
#' @param filenames Vector of filenames to load using readRDS/tbk_data
#' @param databaseSets List of vectors corresponding to probe locations for
#' which the features will be extracted
#'
#' @return DataFrame where rows are sample names and columns are features across
#' different databaseSets
buildStatisticDataSet = function(betas, databaseSets) {
    a = t(data.frame(setNames(
        lapply(filenames,
               function(filename) {
                   unlist(calcDatabaseSetStatisticsAll(betas,
                                                       databaseSets))
               }),
        lapply(filenames,
               function(filename){
                   strsplit(basename(filename), '\\.')[[1]][1]
               }))))
}


#' skew determines the skew of a distribution x, taken from the Moments package
#'
#' @param x Vector of numeric values
#' @param na.rm Logical value corresponding to whether NA will be ignored
#'
#' @return Numeric value quantifying the skew of the distribution x
skew = function (x, na.rm = FALSE) {
    if (is.matrix(x))
        apply(x, 2, skew, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm)
            x <- x[!is.na(x)]
        n <- length(x)
        (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
        vapply(x, skew, na.rm = na.rm)
    else skew(as.vector(x), na.rm = na.rm)
}
