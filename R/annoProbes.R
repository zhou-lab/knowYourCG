#' Annotate Probe IDs using KYCG databases
#'
#' see sesameData_annoProbes if you'd like to annotate by genomic coordinates
#' (in GRanges)
#' @param probeIDs probe IDs in a character vector
#' @param databases character or actual database (i.e. list of probe IDs)
#' @param db_names specific database (default to all databases)
#' @param platform EPIC, MM285 etc. will infer from probe IDs if not given
#' @param indicator return the indicator matrix instead of a concatenated
#' annotation (in the case of have multiple annotations)
#' @param sep delimiter used in paste
#' @param silent suppress message
#' @return named annotation vector, or indicator matrix
#' @examples
#'
#' kycgDataCache(data_titles="KYCG.MM285.designGroup.20210210")
#' probes <- names(sesameData::sesameData_getManifestGRanges("MM285"))
#' anno <- annoProbes(probeIDs=probes, "designGroup", silent = TRUE)
#' @export
annoProbes <- function(probeIDs, databases, db_names = NULL,
    platform = NULL, sep = ",", indicator = FALSE, silent = FALSE) {

    platform <- queryCheckPlatform(platform, probeIDs, silent = silent)
    if (is.character(databases)) {
        dbs <- getDBs(databases, db_names = db_names,
            platform = platform, silent = silent, type = "categorical")
    } else {
        dbs <- databases
    }

    if (is.null(dbs)) { # if there is no dbs, return all-NA
        return(setNames(rep(NA, length(probeIDs)), probeIDs))
    }

    if (is.null(names(dbs)) || any(!nzchar(names(dbs)))) {
        names(dbs) <- paste0("db", seq_along(dbs))
    }

    indicator_mat <- do.call(cbind, lapply(names(dbs), function(db_nm) {
        db <- dbs[[db_nm]]
        probeIDs %in% db
    }))
    if (indicator) {
        rownames(indicator_mat) <- probeIDs
        colnames(indicator_mat) <- names(dbs)
        return(indicator_mat)
    } else {
        anno <- apply(
            indicator_mat, 1,
            function(x) paste(names(dbs)[x], collapse = sep)
        )
        anno <- ifelse(anno == "", NA, anno)
        names(anno) <- probeIDs
        return(anno)
    }
}
