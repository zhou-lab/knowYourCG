#' calcDatabaseSetStatistics1 calculates features of x
#'
#' @param x Vector of numeric values
#'
#' @return Vector with ~20 different engineered features
#'
#' @import stats
calcDatabaseSetStatistics1 = function(x) {
    c(mean=mean(x, na.rm=TRUE),
      median=median(x, na.rm=TRUE),
      var=var(x, na.rm=TRUE),
      stdev=sd(x, na.rm=TRUE),
      skew=as.numeric(skew(x, na.rm=TRUE)),
      iqr=IQR(x, na.rm=TRUE),
      range=max(x, na.rm=TRUE) - min(x, na.rm=TRUE),
      min=min(x, na.rm=TRUE),
      max=max(x, na.rm=TRUE),
      quantile(x, na.rm=TRUE, probs=seq(0, 1, 0.1)))
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
calcDatabaseSetStatisticsAll = function(filename, databaseSets) {
    unlist(lapply(names(databaseSets),
                  function(databaseSetName) {
                      databaseSet = databaseSets[[databaseSetName]]
                      # statistics = calcDatabaseSetStatistics1(tbk_data(filename,
                      # probes=head(databaseSet)))
                      rds = readRDS(filename)
                      statistics = calcDatabaseSetStatistics1(rds[names(databaseSets)])
                      names(statistics) = lapply(names(statistics), function(colname) {
                          paste(databaseSetName, colname, sep="-")
                      })
                      return(statistics)
                  }))
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
buildStatisticDataSet = function(filenames, databaseSets) {
    a = t(data.frame(setNames(
        lapply(filenames,
               function(filename) {
                   unlist(calcDatabaseSetStatisticsAll(filename,
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
