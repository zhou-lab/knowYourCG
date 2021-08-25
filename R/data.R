baseurl = "http://zhouserver.research.chop.edu"

#' databaseSetGet retrieves database sets from a meta data sheet by querying the
#' group, platform, reference columns. The data is returned as a list where the
#' names correspond to chosen databaseSets.
#'
#' @param accessions vector containing the characters associated with the
#' selected databaseSets; only non-NA locations will be returned. Optional.
#' (Default: c("20210810_MM285_TFBS_ENCODE").
#' @param group string representing the group for which the databaseSets will
#' be returned. Optional. (Default: NA).
#' @param platform string representing the platform (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned. Optional. (Default: NA).
#' @param cacheLoc String corresponding to the local filesystem location of
#' where the cache should be stored. Optional. (Default: "").
#' @param release Integer indicating the release number of the databaseSet
#' manifest to use. Optional. (Defualt: 2).
#' @param dev Logical value indiciating whether to use development version
#' of the manifest file. Optional. (Default: TRUE).
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return One list of vectors corresponding to aggregated databaseSets.
#'
#' @examples
#' databaseSetGet()
#'
#' @import readxl
#'
#' @export
databaseSetGet = function(keys=NA, group=NA, platform=NA, reference=NA,
                          cacheLoc="", release=2, dev=TRUE, verbose=TRUE) {
    options(timeout=1000)

    if (verbose) {
        print(sprintf("Loading in databaseSet manifest release %d...", release))
    }

    if (dev) {
        meta = read_excel("/Users/ethanmoyer/Dropbox/Ongoing_knowYourCpG/20210710_databaseSets.xlsx", "R2 In Progress")
        meta = meta[as.logical(meta$Development), ]
    } else {
        meta = read.table(url(sprintf("%s/kyCG/RELEASE_%s.csv",
                                   baseurl, release)), header=TRUE)
    }

    if (any(!is.na(keys))) {
        meta = meta[match(keys, meta$Accession), ]
    }

    if (!is.na(group)) {
        meta = meta[match(group, meta$Group), ]
    }

    if (!is.na(platform)) {
        meta = meta[match(platform, meta$Array), ]
    }

    if (!is.na(reference)) {
        meta = meta[match(reference, meta$Reference), ]
    }

    if (verbose) {
        print(sprintf("Retrieving %d databaseSets...", sum(meta$N)))
    }

    cacheDatabaseSets(release=release)

    path = Sys.getenv("KYCG_DATABASESETS_LOC")
    if (path == "") {
        path = file.path('..', 'databaseSets')
    }

    path = file.path(path, sprintf("RELEASE_%s", release))

    bfc = BiocFileCache(path, ask = FALSE)
    bfcinfoAll = bfcinfo(bfc)

    filelocations = lapply(meta$Accession, function(acc) bfcquery(bfc, acc, field=c("rname"))$rpath)

    databaseSets = flattenlist(lapply(filelocations, function(filelocation) readRDS(filelocation)))

    return(databaseSets)
}


#' flattenlist flattens a multidimensional list into a single dimensional list.
#'
#' @return A single dimensional list.
#'
#' @examples
#' flattenlist(list(a=list(1,2,3), b=list(4,5,6)))
flattenlist = function(x) {
    morelists = vapply(x, function(x_) is(x_, 'list'), TRUE)
    out = c(x[!morelists], unlist(x[morelists], recursive=FALSE))
    if(sum(morelists)){
        Recall(out)
    } else{
        return(out)
    }
}


#' listDatabaseSets prints which databaseSets are available for a given release
#'
#' @param release Integer indicating the release number of the databaseSet
#' manifest to use. Optional. (Defualt: 2).
#' @param dev Logical value indiciating whether to use development version
#' of the manifest file. Optional. (Default: TRUE).
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return One list of vectors corresponding to aggregated databaseSets.
#'
#' @examples
#' databaseSetGet(c(1))
#'
#' @import readxl
#'
#' @export
listDatabaseSets = function(release=2, dev=TRUE, verbose=TRUE) {
    if (dev) {
        meta = read_excel("/Users/ethanmoyer/Dropbox/Ongoing_knowYourCpG/20210710_databaseSets.xlsx",
                          "R2 In Progress")
    } else {
        meta = read.table(url(sprintf("%s/kyCG/RELEASE_%s.csv",
                                      baseurl, release)), header=TRUE)
    }

    x = apply(meta, 1, function(row) {if (dev & !as.logical(row["Development"])) return(NULL)
        cat(sprintf("Accession: %s (n: %s)\n", format(row["Accession"], width = 50, justify = "l"), row["N"]))})
}


#' getQuerySets retrieves sample query sets.
#'
#' @return List of categorical vectors, each of which corresponding to a query
#' set.
#'
#' @examples
#' getQuerySets()
#'
#' @export
getQuerySets = function() {
    return(readRDS(url(sprintf("%s/kyCG/20210726_querySets.rds", baseurl))))
}


#' getUniverseSet retrieves universe set of a given array.
#'
#' @param platform string representing the platform (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return Vector of strings corresponding to the ProbeIDs on a given array.
#'
#' @examples
#' getUniverseSet("MM285.mm10.manifest")
#'
#' @export
getUniverseSet = function(platform, verbose=TRUE) {
    tools::R_user_dir("", which="cache")

    path = Sys.getenv("KYCG_UNIVERSESETS_LOC")
    if (path == "") {
        path = file.path("..", "universeSets")
    }

    bfc = BiocFileCache(path, ask = FALSE)

    ridsAll = bfcinfo(bfc)$rname

    manifests = c("MM285.mm10.manifest",
                  "EPIC.hg19.manifest",
                  "HM450.hg19.manifest",
                  "HM27.hg19.manifest")

    manifest = manifests[grepl(platform, manifests)]

    if (length(manifest) == 0) {
        print(sprintf("Invalid platform [%s]", platform))
        return(NULL)
    }

    if (any(manifest %in% ridsAll))
        return(readRDS(bfcquery(bfc, manifest, field=c("rname"))$rpath))

    universeSet = fread(
        sprintf("%s/InfiniumAnnotation/current/%s/%s.tsv.gz", baseurl, platform, manifest))$probeID

    savepath = bfcnew(bfc, manifest, ext=".RDS")
    saveRDS(universeSet, file=savepath)

    return(universeSet)
}


#' getProbeID2Gene retrieves mapping of of genes to probeID for a specific
#' platform.
#'
#' @return List with each probe mapped to its respective gene.?
#'
#' @examples
#' getProbeID2Gene('EPIC')
getProbeID2Gene = function(platform) {
    # TODO: TRY BLOCK
    tryCatch({
        readRDS(
            url(sprintf("%s/kyCG/20210726_%s_probeID2gene.rds",
                        baseurl, platform)))
    },
    error = function (condition) {
        print("ERROR:")
        print(paste("  Message:",conditionMessage(condition)))
        print(paste("  Call: ",conditionCall(condition)))
    },
    finally= function() {
        print(sprintf("Invalid platform [%s]", platform))
    })
}

#' getGSM2Target retrieves sample query sets.
#'
#' @return List of categorical vectors, each of which corresponding to a query
#' set.
#'
#' @examples
#' getGSM2Target('EPIC', 'TFBS')
# getGSM2Target = function(platform, modality) {
#     rds = readRDS(
#         url(sprintf("%s/kyCG/%s.%s.GSM2Target.rds",
#                     baseurl, platform, modality)))
#     return(rds)
# }



#' cacheDatabaseSets cache databaseSets into memory
#'
#' @param release Integer indicating the release number of the databaseSet
#' manifest to use. Optional. (Defualt: 2).
#' @param dev Logical value indiciating whether to use development version
#' of the manifest file. Optional. (Default: TRUE).
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return One list of vectors corresponding to aggregated databaseSets.
#'
#' @import tools
#' @import dbplyr
#' @import BiocFileCache
#'
#' @examples
#' cacheDatabaseSets(release=2)
cacheDatabaseSets = function(release=2, dev=TRUE, verbose=TRUE) {
    tools::R_user_dir("", which="cache")

    path = Sys.getenv("KYCG_DATABASESETS_LOC")
    if (path == "") {
        path = file.path("..", "databaseSets")
    }

    path = file.path(path, sprintf("RELEASE_%s", release))

    if (dev) {
        meta = read_excel("/Users/ethanmoyer/Dropbox/Ongoing_knowYourCpG/20210710_databaseSets.xlsx",
                          "R2 In Progress")
        meta = meta[as.logical(meta$Development), ]
    } else {
        meta = readRDS(url(sprintf("%s/kyCG/20210710_databaseSets.rds",
                                   baseurl)))
        meta = read_excel(sprintf("%s/kyCG/20210710_databaseSets.xlsx",
                                  baseurl), sprintf("R%d", release))
    }

    remotelocations = unlist(lapply(meta$Accession, function(acc) sprintf("%s/kyCG/%s.rds", baseurl, acc)))

    ridsQuery = setNames(remotelocations, unlist(lapply(
        remotelocations,
        function(location) {
            strsplit(strsplit(location, "/")[[1]][5], ".rds")[[1]]
        })))

    bfc = BiocFileCache(path, ask = FALSE)
    bfcinfoAll = bfcinfo(bfc)

    ridsAll = bfcinfoAll$rname

    rids = ridsQuery[!(names(ridsQuery) %in% ridsAll)]

    nrids = length(rids)

    if (nrids > 0) {
        for (irid in 1:nrids) {
            if (verbose)
                cat(sprintf("Downloading ... [%d of %d].\r", irid, nrids))
            rid = rids[irid]
            location=gsub("https", "http", as.character(rid))
            databaseSet = readRDS(url(location))
            savepath = bfcnew(bfc,
                              names(rid),
                              ext=".RDS")
            saveRDS(databaseSet, file=savepath)
        }
    }
}



#' clearCache removes current cache (if it exists) otherwise returns an error.
#'
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return NULL
#'
#' @import dbplyr
#' @import BiocFileCache
#'
#' @examples
#'
#' @export
clearCache = function(verbose=TRUE) {
    path = Sys.getenv("KYCG_DATABASESETS_LOC")
    if (path == "") {
        path = "../databaseSets"
    }
    bfc = BiocFileCache(path, ask = FALSE)
    if (verbose)
        cat(sprintf("Clearing cache ...\r"))
    removebfc(bfc, ask=FALSE)
    if (verbose)
        cat(sprintf("Cache cleared.\n"))
}


