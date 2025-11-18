#' List database group names
#'
#' @param filter keywords for filtering
#' @param path file path to downloaded knowledgebase sets
#' @param type categorical, numerical (default: all)
#' @return a list of db group names
#' @importFrom stringr str_replace
#' @importFrom utils read.table
#' @examples
#' head(listDBGroups("chromHMM"))
#' ## or listDBGroups(path = "~/Downloads")
#' @export
listDBGroups <- function(filter = NULL, path = NULL, type = NULL) {

    if (is.null(path)) {
        gps <- df_master[,c("Title","Description")]
        gps$type <- vapply(strsplit(
            gps$Description, " "), function(x) x[2], character(1))
        gps$Description <- str_replace(
            gps$Description, "KYCG categorical database holding ", "")
        if (!is.null(filter)) {
            gps <- gps[grepl(filter, gps$Title),]
        }
        if (!is.null(type)) {
            gps <- gps[gps$type %in% type,]
        }
    } else {
        gps <- basename(list.files(path, recursive = TRUE))
    }
    gps
}
