#' Load knowledgebase databases from TSV files
#'
#' This used to be an exported function. Now it's internal. Use RDS download
#' directly.
#'
#' This function loads knowledgebase sets from tab-delimited (.tsv or .tsv.gz)
#' files downloaded from Zenodo or other sources. The TSV files should contain
#' two columns: "Probe_ID" and "Knowledgebase". The function splits the data
#' by knowledgebase name and returns a list of database vectors.
#'
#' @param in_paths Character vector of file paths, URLs to .tsv or
#' .tsv.gz files, or a single directory path containing such files.
#' If a directory is provided, all files in that directory will be loaded.
#' URLs (starting with http:// or https://) are loaded directly.
#' @return A list of database vectors. Each element contains Probe_IDs with
#' attributes:
#' \itemize{
#'   \item \code{group} - The database group name
#'     (derived from filename)
#'   \item \code{dbname} - The knowledgebase name
#'     (from the Knowledgebase column)
#' }
#' @details
#' The input TSV file(s) must have a header row and contain at least
#' two columns:
#' \itemize{
#'   \item \code{Probe_ID} - Probe identifiers (e.g., cg12345678)
#'   \item \code{Knowledgebase} - Knowledgebase/database name
#' }
#' @importFrom readr read_tsv
loadDBs <- function(in_paths) {
    # Detect URLs for naming purposes
    is_url <- grepl("^https?://", in_paths)

    if (length(in_paths)==1 && !any(is_url) && dir.exists(in_paths)) {
        groupnms <- list.files(in_paths)
        in_paths <- file.path(in_paths, groupnms)
    } else {
        groupnms <- basename(in_paths)
    }

    dbs <- do.call(c, lapply(seq_along(groupnms), function(i) {
        # readr::read_tsv handles both URLs and local files
        tbl <- readr::read_tsv(in_paths[i], show_col_types = FALSE)
        dbs1 <- split(tbl$Probe_ID, tbl$Knowledgebase)
        lapply(names(dbs1), function(dbname) {
            db1 <- dbs1[[dbname]];
            attr(db1, "group") <- sub(".gz$","",groupnms[i]);
            attr(db1, "dbname") <- dbname;
            db1;})
    }))

    names(dbs) <- vapply(dbs, function(x) attr(x, "dbname"), character(1))
    dbs
}

