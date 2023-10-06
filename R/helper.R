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

inferUniverse <- function(platform) {
    mfts <- c(
        "MM285.address", "EPIC.address",
        "Mammal40.address", "HM450.address", "HM27.address")
    mft <- mfts[grepl(platform, mfts)]
    stopifnot(length(mft) == 1 && all(mft %in% mfts))
    sesameDataGet(mft)$ordering$Probe_ID
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
        nms, platform = NULL,allow_multi = FALSE, type = NULL,
        silent = FALSE) {

    gps <- KYCG_listDBGroups(type = type)
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
        message("Selected the following database groups:")
        invisible(lapply(seq_along(nms), function(i) {
            message(sprintf("%d. %s", i, nms[i]))
        }))
    }
    nms
}


#' List database group names
#'
#' @param filter keywords for filtering
#' @param path file path to downloaded knowledgebase sets
#' @param type categorical, numerical (default: all)
#' @return a list of db group names
#' @examples
#' head(KYCG_listDBGroups("chromHMM"))
#' ## or KYCG_listDBGroups(path = "~/Downloads")
#' @export
KYCG_listDBGroups <- function(filter = NULL, path = NULL, type = NULL) {

    if (is.null(path)) {
        gps <- sesameDataList("KYCG", full=TRUE)[,c("Title","Description")]
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



KYCG_loadDBs <- function(in_paths) {
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



#' Get databases by full or partial names of the database group(s)
#'
#' @param group_nms database group names
#' @param db_names name of the database, fetech only the given databases
#' @param platform EPIC, HM450, MM285, ... If given, will restrict to
#' that platform.
#' @param summary return a summary of database instead of db itself
#' @param allow_multi allow multiple groups to be returned for
#' @param type numerical, categorical, default: all
#' @param silent no messages
#' each query.
#' @return a list of databases, return NULL if no database is found
#' @examples
#' dbs <- KYCG_getDBs("MM285.chromHMM")
#' dbs <- KYCG_getDBs(c("MM285.chromHMM", "MM285.probeType"))
#' @export
KYCG_getDBs <- function(
        group_nms, db_names = NULL, platform = NULL,
        summary = FALSE, allow_multi = FALSE, type = NULL, silent = FALSE) {

    if (!is.character(group_nms)) {
        return(group_nms)
    }

    group_nms <- guess_dbnames(group_nms, platform = platform,
                               allow_multi = TRUE, type = type,
                               silent = silent)
    ## group_nms <- group_nms[sesameDataHas(group_nms)]
    group_nms <- group_nms[group_nms %in% sesameDataList()$Title]
    if (length(group_nms) == 0) {
        return(NULL)
    }
    res <- do.call(c, lapply(unname(group_nms), function(nm) {
        dbs <- sesameDataGet(nm)
        setNames(lapply(seq_along(dbs), function(ii) {
            db <- dbs[[ii]]
            attr(db, "group") <- nm
            attr(db, "dbname") <- names(dbs)[ii]
            db
        }), names(dbs))}))

    if (summary) {
        do.call(bind_rows, lapply(res, attributes))
    } else if (is.null(db_names)) {
        res
    } else {
        stopifnot(all(db_names %in% names(res)))
        res[db_names]
    }
}



#' Annotate Probe IDs using KYCG databases
#'
#' see sesameData_annoProbes if you'd like to annotate by genomic coordinates
#' (in GRanges)
#' @param query probe IDs in a character vector
#' @param databases character or actual database (i.e. list of probe IDs)
#' @param db_names specific database (default to all databases)
#' @param platform EPIC, MM285 etc. will infer from probe IDs if not given
#' @param indicator return the indicator matrix instead of a concatenated
#' annotation (in the case of have multiple annotations)
#' @param sep delimiter used in paste
#' @param silent suppress message
#' @return named annotation vector, or indicator matrix
#' @examples
#' query <- names(sesameData_getManifestGRanges("MM285"))
#' anno <- KYCG_annoProbes(query, "designGroup", silent = TRUE)
#' @export
KYCG_annoProbes <- function(query, databases, db_names = NULL,
                            platform = NULL, sep = ",",
                            indicator = FALSE, silent = FALSE) {

    platform <- queryCheckPlatform(platform, query, silent = silent)
    if (is.character(databases)) {
        dbs <- KYCG_getDBs(databases, db_names = db_names,
                           platform = platform, silent = silent,
                           type = "categorical")
    } else {
        dbs <- databases
    }

    ind <- do.call(cbind, lapply(names(dbs), function(db_nm) {
        db <- dbs[[db_nm]]
        query %in% db
    }))
    if (indicator) {
        rownames(ind) <- query
        colnames(ind) <- names(dbs)
        return(ind)
    } else {
        anno <- apply(ind, 1, function(x) paste(names(dbs)[x], collapse=sep))
        anno <- ifelse(anno == "", NA, anno)
        names(anno) <- query
        return(anno)
    }
}
