#' Load knowledgebase databases from TSV files
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
#' @examples
#' 
#' # Load directly from a URL
#' dbs <- loadDBs(
#'   "https://zenodo.org/records/18176501/files/ImprintingDMR.20220818.gz")
#'
#' # Examine the structure
#' length(dbs)  # Number of databases loaded
#' names(dbs)   # Database names
#' head(dbs[[1]])  # First database content
#'
#' # Load from multiple URLs
#' urls <- c(
#'   "https://zenodo.org/records/18176501/files/ImprintingDMR.20220818.gz",
#'   "https://zenodo.org/records/18176501/files/Blacklist.20220304.gz"
#' )
#' dbs_multi <- loadDBs(urls)
#'
#' # Load from local file (download to temp file first)
#' tmp_file <- tempfile(fileext = ".tsv.gz")
#' download.file(
#'   "https://zenodo.org/records/18176501/files/ImprintingDMR.20220818.gz",
#'   tmp_file, mode = "wb", quiet = TRUE
#' )
#' dbs_local <- loadDBs(tmp_file)
#' unlink(tmp_file)  # Clean up
#'
#' # Load all files from a directory
#' tmp_dir <- tempfile()
#' dir.create(tmp_dir)
#' download.file(
#'   "https://zenodo.org/records/18176501/files/ImprintingDMR.20220818.gz",
#'   file.path(tmp_dir, "ImprintingDMR.20220818.gz"), mode = "wb", quiet = TRUE
#' )
#' download.file(
#'   "https://zenodo.org/records/18176501/files/Blacklist.20220304.gz",
#'   file.path(tmp_dir, "Blacklist.20220304.gz"), mode = "wb", quiet = TRUE
#' )
#' dbs_all <- loadDBs(tmp_dir)
#' unlink(tmp_dir, recursive = TRUE)  # Clean up
#'
#' @importFrom readr read_tsv
#' @export
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

