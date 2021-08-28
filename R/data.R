baseurl = "http://zhouserver.research.chop.edu"

#' getDatabaseSets retrieves database sets from a meta data sheet by querying the
#' group, array, reference columns. The data is returned as a list where the
#' names correspond to chosen databaseSets.
#'
#' @param accessions vector containing the characters associated with the
#' selected databaseSets; only non-NA locations will be returned. Optional.
#' (Default: c("20210810_MM285_TFBS_ENCODE").
#' @param group string representing the group for which the databaseSets will
#' be returned. Optional. (Default: NA).
#' @param array string representing the array (EPIC, HM450, HM27, MM285)
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
#' getDatabaseSets()
#'
#' @import readxl
#'
#' @export
getDatabaseSets = function(keys=NA, group=NA, array=NA, reference=NA,
                          cacheLoc="", release=2, dev=TRUE, verbose=TRUE) {
    options(timeout=1000)

    if (verbose) {
        print(sprintf("Loading in databaseSet manifest release %d...", release))
    }

    if (dev) {
        meta = read_excel(sprintf("%s/Dropbox/Ongoing_knowYourCpG/20210710_databaseSets.xlsx", Sys.getenv("HOME")), "R2 In Progress")
        meta = meta[as.logical(meta$Development), ]
    } else {
        meta = read.table(url(sprintf("%s/kyCG/RELEASE_%s.csv",
                                   baseurl, release)), header=TRUE)
    }

    if (any(!is.na(keys))) {
        meta = meta[which(meta$Accession %in% keys), ]
    }
    
    if (!is.na(group)) {
        meta = meta[which(meta$Group %in% group), ]
    }

    if (!is.na(array)) {
        meta = meta[which(meta$Array %in% array), ]
    }

    if (!is.na(reference)) {
        meta = meta[which(meta$Reference %in% reference), ]
    }
    
    if (verbose) {
        print(sprintf("Retrieving %d databaseSets...", sum(meta$N)))
    }

    cacheDatabaseSets(release=release)

    path = Sys.getenv("KYCG_DATABASESETS_LOC")
    if (path == "") {
        path = file.path('databaseSets')
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
#' listDatabaseSets(release=2)
#'
#' @import readxl
#'
#' @export
listDatabaseSets = function(release=2, dev=TRUE, verbose=TRUE) {
    if (dev) {
        meta = read_excel(sprintf("%s/Dropbox/Ongoing_knowYourCpG/20210710_databaseSets.xlsx", Sys.getenv("HOME")),
                          "R2 In Progress", )
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
#' @param array string representing the array (EPIC, HM450, HM27, MM285)
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
getUniverseSet = function(array, verbose=TRUE) {
    tools::R_user_dir("", which="cache")

    path = Sys.getenv("KYCG_UNIVERSESETS_LOC")
    if (path == "") {
        path = file.path("universeSets")
    }

    bfc = BiocFileCache(path, ask = FALSE)

    ridsAll = bfcinfo(bfc)$rname

    manifests = c("MM285.mm10.manifest",
                  "EPIC.hg19.manifest",
                  "HM450.hg19.manifest",
                  "HM27.hg19.manifest")

    manifest = manifests[grepl(array, manifests)]

    if (length(manifest) == 0) {
        print(sprintf("Invalid array [%s]", array))
        return(NULL)
    }
    
    if (length(grepl(manifest, ridsAll)) == 1 && grepl(manifest, ridsAll))
        return(readRDS(bfcquery(bfc, manifest, field=c("rname"))$rpath))

    manifest_data = fread(
        sprintf("%s/InfiniumAnnotation/current/%s/%s.tsv.gz", baseurl, array, manifest))
    
    universeSet = manifest_data$probeID
    
    if (is.null(universeSet))
        universeSet = manifest_data$Probe_ID

    savepath = bfcnew(bfc, manifest, ext=".RDS")
    saveRDS(universeSet, file=savepath)

    return(universeSet)
}


#' getProbeID2Gene retrieves mapping of of genes to probeID for a specific
#' array.
#' 
#' @param array string representing the array (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#'
#' @return List with each probe mapped to its respective gene.?
#'
#' @examples
#' getProbeID2Gene('EPIC')
getProbeID2Gene = function(array) {
    # TODO: TRY BLOCK
    tryCatch({
        readRDS(
            url(sprintf("%s/kyCG/20210726_%s_probeID2gene.rds",
                        baseurl, array)))
    },
    error = function (condition) {
        print("ERROR:")
        print(paste("  Message:",conditionMessage(condition)))
        print(paste("  Call: ",conditionCall(condition)))
    },
    finally= function() {
        print(sprintf("Invalid array [%s]", array))
    })
}


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
        path = file.path("databaseSets")
    }

    path = file.path(path, sprintf("RELEASE_%s", release))

    if (dev) {
        meta = read_excel(sprintf("%s/Dropbox/Ongoing_knowYourCpG/20210710_databaseSets.xlsx", Sys.getenv("HOME")),
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



#' clearCache removes all current cache (if it exists).
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
#' clearCache()
#'
#' @export
clearCache = function(verbose=TRUE) {
    clearDatabaseSetCache(verbose=verbose)
    clearUniverseSetCache(verbose=verbose)
    clearFeatureEngineeringCache(verbose=verbose)
}

#' clearDatabaseSetCache removes database set current cache (if it exists).
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
#' clearDatabaseSetCache()
#'
#' @export
clearDatabaseSetCache = function(verbose=TRUE) {
    if (verbose)
        cat(sprintf("Clearing database set cache ...\r"))
    path = Sys.getenv("KYCG_DATABASESETS_LOC")
    if (path == "") {
        path = "databaseSets"
    }
    
    bfc = BiocFileCache(path, ask = FALSE)
    removebfc(bfc, ask=FALSE)
    
    if (verbose)
        cat(sprintf("Cache database set cleared.\n"))
}

#' clearUniverseSetCache removes universe set current cache (if it exists).
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
#' clearUniverseSetCache()
#'
#' @export
clearUniverseSetCache = function(verbose=TRUE) {
    if (verbose)
        cat(sprintf("Clearing universe set cache ...\r"))
    path = Sys.getenv("KYCG_UNIVERSESETS_LOC")
    if (path == "") {
        path = "universeSets"
    }
    bfc = BiocFileCache(path, ask = FALSE)
    removebfc(bfc, ask=FALSE)

    if (verbose)
        cat(sprintf("Cache universe set cleared.\n"))
}

#' clearFeatureEngineeringCache removes universe set current cache (if it 
#' exists).
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
#' clearFeatureEngineeringCache()
#'
#' @export
clearFeatureEngineeringCache = function(verbose=TRUE) {
    if (verbose)
        cat(sprintf("Clearing feature engineering cache ...\r"))
    
    path = Sys.getenv("KYCG_FEATURE_ENGINEERING_LOC")
    if (path == "") {
        path = "featureEngineering"
    }
    bfc = BiocFileCache(path, ask = FALSE)
    removebfc(bfc, ask=FALSE)
    
    if (verbose)
        cat(sprintf("Cache feature engineering cleared.\n"))
}


#' getSampleSheet retrieves an example samplesheet composed of 
#' 20 samples.
#'
#' @param array string representing the array (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#' 
#' @param array string representing the array (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#' @param release Integer indicating the release number of the databaseSet
#' manifest to use. Optional. (Defualt: 2).
#' @param dev Logical value indiciating whether to use development version
#' of the manifest file. Optional. (Default: TRUE).
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return samplesheet
#'
#' @import dbplyr
#' @import BiocFileCache
#'
#' @examples
#' getSampleSheet("MM285")
#'
#' @export
getSampleSheet = function(array, release=1, dev=TRUE, verbose=TRUE) {
    if (verbose) {
        print(sprintf("Loading in Mouse samplesheet release %d...", release))
    }
    
    if (dev) {
        meta = read_excel('/Users/moyerej/Dropbox/Ongoing_knowYourCpG/20210827_featureEngineering.xlsx', sprintf("R%s", release))
    } else {
        meta = NULL
    }
        
    path = Sys.getenv("KYCG_FEATURE_ENGINEERING_LOC")
    if (path == "") {
        path = file.path("featureEngineering")
    }
    
    bfc = BiocFileCache(path, ask = FALSE)
    
    ridsAll = bfcinfo(bfc)$rname
    
    query = sprintf("%s_samplesheet", toupper(array))
    
    if (length(grepl(query, ridsAll)) == 1 && grepl(query, ridsAll))
        return(readRDS(bfcquery(bfc, query, field=c("rname"))$rpath))
    
    meta = meta[which(meta$Group %in% "samplesheet"), ]

    meta = meta[which(meta$Array %in% array), ]
    
    samplesheet = readRDS(url(sprintf("%s/kyCG/%s", baseurl, meta$Accession)))
    
    savepath = bfcnew(bfc, meta$Accession, ext=".RDS")
    saveRDS(samplesheet, file=savepath)
    
    return(samplesheet)
}

#' getBetas retrieves an example betas matrix composed of 20 samples.
#'
#' @param array string representing the array (EPIC, HM450, HM27, MM285)
#' for which databaseSets will be returned.
#' @param release Integer indicating the release number of the databaseSet
#' manifest to use. Optional. (Defualt: 2).
#' @param dev Logical value indiciating whether to use development version
#' of the manifest file. Optional. (Default: TRUE).
#' @param verbose Logical value indicating whether intermediate outputs will be
#' displayed to console. Optional. (Default: TRUE).
#'
#' @return betas matrix
#'
#' @import dbplyr
#' @import BiocFileCache
#'
#' @examples
#' getBetas("MM285")
#'
#' @export
getBetas = function(array, release=1, dev=TRUE, verbose=TRUE) {
    if (verbose) {
        print(sprintf("Loading in Mouse samplesheet release %d...", release))
    }
    
    if (dev) {
        meta = read_excel('/Users/moyerej/Dropbox/Ongoing_knowYourCpG/20210827_featureEngineering.xlsx', sprintf("R%s", release))
    } else {
        meta = NULL
    }
    
    path = Sys.getenv("KYCG_FEATURE_ENGINEERING_LOC")
    if (path == "") {
        path = file.path("featureEngineering")
    }
    
    bfc = BiocFileCache(path, ask = FALSE)
    
    ridsAll = bfcinfo(bfc)$rname
    
    query = sprintf("%s_betas", toupper(array))
    
    if (length(grepl(query, ridsAll)) == 1 && grepl(query, ridsAll))
        return(readRDS(bfcquery(bfc, query, field=c("rname"))$rpath))
    
    meta = meta[which(meta$Group %in% "betas"), ]
    
    meta = meta[which(meta$Array %in% array), ]
    
    betas = readRDS(url(sprintf("%s/kyCG/%s", baseurl, meta$Accession)))
    
    savepath = bfcnew(bfc, meta$Accession, ext=".RDS")
    saveRDS(betas, file=savepath)
    
    return(betas)
}

