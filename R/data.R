baseurl = "http://zhouserver.research.chop.edu"

#' databaseSetGet retrieves database sets from a meta data sheet by querying the array, group, and reference columns. The data is returned as a list where the names correspond to chosen database sets.
#'
#' @param keys vector containing the keys associated with the selected
#' databaseSets; only non-NA locations will be returned.
#' @param group string representing the group for which the databaseSets will
#' be returned.
#' @param platform string representing the platform (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#'
#' @return One list of vectors corresponding to aggregated databaseSets.
#'
#' @examples
#' databaseSetGet(c(1))
#'
#' @export
databaseSetGet = function(keys=-1, group=NA, platform=NA, cacheLoc="") {
    options(timeout=1000)
    meta = readRDS(url(sprintf("%s/kyCG/20210710_databaseSets.rds", baseurl)))

    if (keys >= 0) {
        meta = meta[match(keys, meta$Key), ]
    }

    if (!is.na(group)) {
        meta = meta[match(group, meta$Group), ]
    }

    if (!is.na(platform)) {
        meta = meta[match(platform, meta$Array), ]
    }

    locations = na.omit(meta$Location)

    if (length(locations) == 0) {
        return(NULL)
    }

    databaseSets = cacheDatabaseSets(locations, cacheLoc=cacheLoc)

    return(databaseSets)
}


#' flattenlist flattens a multidimensional list into a single dimensional list.
#'
#' @return A single dimensional list.
#'
#' @examples
#' flattenlist(list(a=list(1,2,3), b=list(4,5,6)))
flattenlist = function(x){
    morelists = vapply(x, function(x_) is(x_, 'list'), TRUE)
    out = c(x[!morelists], unlist(x[morelists], recursive=FALSE))
    if(sum(morelists)){
        Recall(out)
     } else{
        return(out)
    }
}

#' getTestQuerySet retrieves sample query sets.
#'
#' @return List of categorical vectors, each of which corresponding to a query
#' set.
#'
#' @examples
#' getTestQuerySet()
#'
#' @export
getTestQuerySet = function() {
    return(readRDS(url(sprintf("%s/kyCG/20210726_testResults.rds", baseurl))))
}

#' getUniverseSet retrieves universe set of a given array.
#'
#' @param platform string representing the platform (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#'
#' @return Vector of strings corresponding to the ProbeIDs on a given array.
#'
#' @examples
#' getUniverseSet("MM285")
#'
#' @export
getUniverseSet = function(platform, cacheLoc="") {
    tools::R_user_dir(cacheLoc, which="cache")

    path = Sys.getenv("kyCpGCache")
    if (path == "") {
        path = tempfile()
        Sys.setenv(kyCpGCache=path)
    }

    bfc = BiocFileCache(path, ask = FALSE)

    ridsAll = bfcinfo(bfc)$rname

    if (any(platform %in% ridsAll))
        return(readRDS(bfcquery(bfc, platform, exact=TRUE)$rpath))

    if (platform == "MM285") {
        universeSet = get(load(url(sprintf("%s/sesameData/MM285.mm10.manifest.rda",
                                           baseurl))))$Probe_ID
    } else if (platform == "EPIC") {
        universeSet = readRDS(url(sprintf("%s/InfiniumAnnotation/current/EPIC/EPIC.hg19.manifest.rds",
                                          baseurl)))$probeID
    } else if (platform == "HM450") {
        universeSet = readRDS(url(sprintf("%s/InfiniumAnnotation/current/HM450/HM450.hg19.manifest.rds",
                                          baseurl)))$probeID
    } else if (platform == "HM27") {
        universeSet = readRDS(url(sprintf("%s/InfiniumAnnotation/current/HM27/HM27.hg19.manifest.rds",
                                          baseurl)))$probeID
    } else {
        print(sprintf("Invalid platform [%s]", platform))
        return(NULL)
    }

    # TODO: Change ID to reflect file name instead of platform

    savepath = bfcnew(bfc, platform, ext=".RDS")
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
        readRDS(url(sprintf("%s/kyCG/20210726_%s_probeID2gene.rds", baseurl, platform)))
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
getGSM2Target = function(platform, modality) {
    return(readRDS(url(sprintf("%s/kyCG/%s.%s.GSM2Target.rds", baseurl, platform, modality))))
}

#' cacheDatabaseSets cache databaseSets into memory
#'
#' @param locations Vector of strings corresponding to the URL locations of
#' databaseSets stored in RDS format.
#' @param cacheLoc String corresponding to the local filesystem location of
#' where the cache should be stored.
#'
#' @return One list of vectors corresponding to aggregated databaseSets.
#'
#' @import tools
#' @import dbplyr
#' @import BiocFileCache
#'
#' @examples
#' cacheDatabaseSets(c("https://zhouserver.research.chop.edu/kyCG/20210713_EPIC_hg19_SNPs.rds", "https://zhouserver.research.chop.edu/kyCG/20210630_MM285_mm10_CTCF.rds"))
cacheDatabaseSets = function(locations, cacheLoc="") {
    tools::R_user_dir(cacheLoc, which="cache")

    path = Sys.getenv("kyCpGCache")
    if (path == "") {
        path = tempfile()
        Sys.setenv(kyCpGCache=path)
    }

    bfc = BiocFileCache(path, ask = FALSE)

    ridsQuery = setNames(locations, unlist(lapply(
        locations,
        function(location) {
            strsplit(location, "/")[[1]][5]
        })))

    ridsAll = bfcinfo(bfc)$rname

    rids = ridsQuery[!(names(ridsQuery) %in% ridsAll)]

    if (length(rids) > 0)
        databaseSetsNew = flattenlist(lapply(as.vector(rids),
                                      function(location) {
                                          location=gsub("https", "http", location)
                                          databaseSet = readRDS(url(location))
                                          savepath = bfcnew(bfc, strsplit(location, "/")[[1]][5], ext=".RDS")
                                          saveRDS(databaseSet, file=savepath)
                                          return(databaseSet)
                                          }))
    else
        databaseSetsNew = c()

    rids = ridsQuery[names(ridsQuery) %in% ridsAll]

    databaseSetsCached = flattenlist(lapply(names(rids),
                                     function(rname) {
                                         readRDS(bfcquery(bfc, rname)$rpath)
                                     }))

    databaseSets = append(databaseSetsNew, databaseSetsCached)

    return(databaseSets)

}
