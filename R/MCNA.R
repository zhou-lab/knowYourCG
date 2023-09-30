imputeRowMean <- function(mtx) {
    k <- which(is.na(mtx), arr.ind=TRUE)
    mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
    mtx
}


cleanMatrix <- function(mtx, f_row = 0.5, f_col = 0.5) {
    message("Before: ", nrow(mtx), " rows\n")
    namtx <- !is.na(mtx)
    good_row <- rowSums(namtx) > 1
    message("After: ", sum(good_row), " rows\n")
    mtx <- mtx[good_row,]
    imputeRowMean(mtx)
}


returnDiffCpGs <- function(
        betas, query, k=50, metric="correlation", diffThreshold=0.5) {
    refGraph <- rnndescent::nnd_knn(betas, k = k, metric=metric)
    searchGraph <- rnndescent::prepare_search_graph(
        betas,
        refGraph,
        metric=metric,
        verbose = FALSE
    )
    query_nn <- rnndescent::graph_knn_query(
        query = query,
        reference = betas,
        reference_graph = searchGraph,
        k = 1,
        metric = metric,
        verbose = FALSE
    )
    query[which(query_nn$dist > diffThreshold), ]
}


prepareSampleSet <- function(
        betas,k=50,impute=TRUE,num_row=75000,diffThreshold=0.5) {

    if (is(betas, "numeric")) {
        betas <- cbind(sample = betas)
    }
    if (impute) {
        betas <- cleanMatrix(betas)
    }
    num <- ifelse(nrow(betas) < num_row,nrow(betas),num_row)
    var_rows <- order(-apply(betas,1,sd))[seq(num)]
    betas <- betas[var_rows,]
    sample_size <- round(.33 * num)
    betas_sample <- betas[sample(rownames(betas), size=sample_size), ]
    query <- betas[!rownames(betas) %in% rownames(betas_sample),]
    betas_sample <- rbind(
        betas_sample,
        returnDiffCpGs(
            betas=betas_sample,
            query=query,
            k=k,
            diffThreshold = diffThreshold
            )
        )
    betas_sample
}


detectCommunity <- function(el,edgeThreshold=.1,nodeThreshold=0) {
    g <- igraph::graph_from_data_frame(
        el,
        directed = FALSE
    )
    g <- igraph::delete.edges(
        g,
        which(igraph::E(g)$dist > edgeThreshold)
    )
    isolated <- which(igraph::degree(g)==nodeThreshold)
    g <- igraph::delete.vertices(g, isolated)
    lc <- cluster_louvain(g)
    lc
}


#' findCpGModules identifies modules of co-methylated CpGs
#'
#' @param betas matrix of beta values where probes are on the rows and
#' samples on the columns
#' @param k # of neighbors to return from reference graph for query CpGs
#' @param diffThreshold Distance to nearest neighbor to determine if
#' query gets added to reference graph
#' @param impute whether to impute missing values using the row mean
#' (Default: TRUE)
#' @param edgeThreshold minimum inter - CpG distance threshold for community
#' detection (1 - correlation)
#' @param nodeThreshold minimum node degree for removal from graph
#' @param metric metric for computing neighbor distance (Default: correlation)
#' @param moduleSize minimum number of CpGs for module consideration
#' @return A list of CpG modules
#' @importFrom igraph graph_from_data_frame delete.edges delete.vertices
#' @importFrom igraph cluster_louvain degree communities sizes
#' @examples
#' library(SummarizedExperiment)
#' se <- sesameDataGet('MM285.467.SE.tissue20Kprobes')
#' betas <- assay(se)
#' modules <- findCpGModules(betas)
#' sesameDataGet_resetEnv()
#'
#' @export
findCpGModules <- function (
        betas,impute=TRUE,diffThreshold=.5,k=50,metric="correlation",
        edgeThreshold=.1,nodeThreshold=0,moduleSize = 5) {

    beta_sample <- prepareSampleSet(
        betas=betas,
        impute=impute,
        diffThreshold=diffThreshold
    )
    nnr <- rnndescent::nnd_knn(beta_sample, k = k, metric=metric)
    nbr_mtx <- nnr$idx[,-1]
    nbrs <- as.vector(nbr_mtx)
    dist <- as.vector(nnr$dist[,-1])
    el <- matrix(0, nrow = nrow(nbr_mtx) * ncol(nbr_mtx), ncol = 2)
    el[,2] <- nbrs
    el[,1] <- rep(seq(nrow(nbr_mtx)), times=ncol(nbr_mtx))
    el_df <- as.data.frame(el);
    select <- !duplicated(t(apply(el_df,1,sort)))
    el_df <- el_df[select,]
    el_df$dist <- dist[select]
    lc <- detectCommunity(
        el=el_df,
        edgeThreshold = edgeThreshold,
        nodeThreshold = nodeThreshold
    )
    communities <- igraph::communities(lc)[igraph::sizes(lc) >= moduleSize]
    modules <- lapply(communities,function(x) {
        indices <- as.numeric(x)
        rownames(beta_sample)[indices]
    })
    names(modules) <- as.character(seq(modules))
    modules
}
