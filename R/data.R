#' databaseSetGet retrieves database sets from a meta data sheet by querying the array, group, and reference columns. The data is returned as a list where the names correspond to chosen database sets.
#'
#' @param keys vector containing the keys associated with the selected
#' databaseSets; only non-NA locations will be returned.
#' @param platform string representing the platform (EPIC, HM450, HM27, MM285)
#' for which database sets will be returned.
#'
#' @return One list containing features corresponding the test estimate, p-value, and type of test.
#'
#' @examples
#' databaseSetGet(c(1))
#'
#' @export
databaseSetGet = function(keys=-1, platform=NA) {
    options(timeout=1000)
    meta = readRDS(url("http://zhouserver.research.chop.edu/
                       kyCG/20210710_databaseSets.rds"))

    if (keys >= 0) {
        meta = meta[match(keys, meta$Key), ]
    }

    if (!is.na(platform)) {
        locations = na.omit(meta$Location[meta$Array == platform])
    }

    databaseSets = flattenlist(lapply(locations,
                                      function(x) {
                                          x=gsub("https", "http", x);
                                          readRDS(url(x))}
                                      ))

    return(databaseSets)
}

flattenlist <- function(x){
    morelists <- vapply(x, function(xprime) class(xprime)[1]=="list")
    out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
    if(sum(morelists)){
        Recall(out)
     } else{
        return(out)
    }
}
