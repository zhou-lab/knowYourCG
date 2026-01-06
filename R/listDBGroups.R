#' List database group names
#'
#' @param filter keywords for filtering
#' @param path file path to downloaded knowledgebase sets
#' @return a list of db group names
#' @importFrom stringr str_replace
#' @importFrom utils read.table
#' @examples
#' head(listDBGroups("chromHMM"))
#' ## or listDBGroups(path = "~/Downloads")
#' @export
listDBGroups <- function(filter = NULL, path = NULL) {

    if (is.null(path)) {
        gps <- df_master[,c("Title","Description")]
        gps$Description <- str_replace(
            gps$Description, "KYCG categorical database holding ", "")
        gps$Description <- str_replace(
            gps$Description, "KYCG knowledgebase holding ", "")
        if (!is.null(filter)) {
            gps <- gps[grepl(filter, gps$Title),]
        }
    } else {
        gps <- basename(list.files(path, recursive = TRUE))
    }
    gps
}
