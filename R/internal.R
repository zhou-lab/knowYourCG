databases_getMeta <- function(dbs) {
    meta <- do.call(bind_rows, lapply(dbs, function(db) {
        m1 <- attributes(db)
        m1 <- m1[names(m1) != "names"] # a special attribute for continuous db
        if ("meta" %in% names(m1)) { # backward compatibility, to delete
            m1 <- c(m1[!(names(m1) %in% c("meta"))], m1$meta)
        } else {
            m1 <- m1[!(names(m1) %in% c("meta"))]
        }

        if (is.null(m1)) {
            data.frame(hasMeta = FALSE)
        } else {
            c(m1, hasMeta = TRUE)
        }
    }))
    meta[,colnames(meta)!="hasMeta"]
}

queryCheckPlatform <- function(platform, query = NULL, silent = FALSE) {
    if (is.null(platform)) {
        stopifnot(!is.null(query))
        if (is.numeric(query)) {
            platform <- inferPlatformFromProbeIDs(
                names(query),
                silent = silent
            )
        } else {
            platform <- inferPlatformFromProbeIDs(query, silent = silent)
        }
    }
    platform
}

subsetDBs <- function(dbs, universe) {
    dbs <- lapply(dbs, function(db) {
        db1 <- intersect(db, universe)
        attributes(db1) <- attributes(db)
        db1
    })
    dbs <- dbs[length(dbs) > 0]
}


guess_dbnames <- function(
        nms, platform = NULL,allow_multi = FALSE, silent = FALSE) {

    gps <- listDBGroups()
    nms <- do.call(c, lapply(nms, function(nm) {
        if (nm %in% gps$Title) {
            return(nm)
        } else if (length(grep(nm, gps$Title)) >= 1) {
            ret <- grep(nm, gps$Title, value=TRUE)
            if (!allow_multi) { ret <- ret[1]; }
            return(ret)
        } else if (length(grep(nm, gps$Title)) == 0) {
            res <- gps$Title[apply(do.call(cbind, lapply(
                strsplit(nm, "\\.")[[1]], function(q1) grepl(q1, gps$Title))),
                1, all)]
            if (length(res) == 1) {
                return(res[1])
            }
        }
        return(nm)
    }))
    if (!is.null(platform)) {
        nms <- grep(platform, nms, value = TRUE)
    }
    if (!silent) {
        if (length(nms) == 0) {
            message("No knowledgebase selected. Please reselect.")
        } else {
            message("Selected the following database groups:")
            invisible(lapply(seq_along(nms), function(i) {
                message(sprintf("%d. %s", i, nms[i]))
            }))
        }
    }
    nms
}

#' Load knowledgebase databases from TSV files
#'
#' This function loads knowledgebase sets from tab-delimited (.tsv or .tsv.gz)
#' files downloaded from Zenodo or other sources. The TSV files should contain
#' two columns: "Probe_ID" and "Knowledgebase". The function splits the data
#' by knowledgebase name and returns a list of database vectors.
#'
#' @param in_paths Character vector of file paths to .tsv or .tsv.gz files,
#' or a single directory path containing such files. If a directory is provided,
#' all files in that directory will be loaded.
#' @return A list of database vectors. Each element contains Probe_IDs with
#' attributes:
#' \itemize{
#'   \item \code{group} - The database group name (derived from filename)
#'   \item \code{dbname} - The knowledgebase name (from the Knowledgebase column)
#' }
#' @details
#' The input TSV file(s) must have a header row and contain at least two columns:
#' \itemize{
#'   \item \code{Probe_ID} - Probe identifiers (e.g., cg12345678)
#'   \item \code{Knowledgebase} - Knowledgebase/database name
#' }
#'
#' @examples
#' \donttest{
#' # Download a knowledgebase TSV file from Zenodo
#' temp_dir <- tempdir()
#' tsv_file <- file.path(temp_dir, "ChromHMM.20220303.gz")
#' download.file(
#'   "https://zenodo.org/records/18176501/files/ChromHMM.20220303.gz",
#'   destfile = tsv_file
#' )
#'
#' # Load the databases from the TSV file
#' dbs <- loadDBs(tsv_file)
#'
#' # Examine the structure
#' length(dbs)  # Number of databases loaded
#' names(dbs)   # Database names
#' head(dbs[[1]])  # First database content
#'
#' # Load multiple files from a directory
#' dbs_all <- loadDBs(temp_dir)
#' }
#'
#' @export
loadDBs <- function(in_paths) {
    if (length(in_paths)==1 && dir.exists(in_paths)) {
        groupnms <- list.files(in_paths)
        in_paths <- file.path(in_paths, groupnms)
    } else {
        groupnms <- basename(in_paths)
    }
    do.call(c, lapply(seq_along(groupnms), function(i) {
        tbl <- read.table(in_paths[i],sep="\t",header=TRUE)
        dbs <- split(tbl$Probe_ID, tbl$Knowledgebase)
        lapply(names(dbs), function(dbname) {
            db1 <- dbs[[dbname]];
            attr(db1, "group") <- sub(".gz$","",groupnms[i]);
            attr(db1, "dbname") <- dbname;
            db1;})
    }))
}

