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

## for loading from .tsv.gz files, internal use
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

