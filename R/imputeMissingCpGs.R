imputeRowMean <- function(mtx) {
    k <- which(is.na(mtx), arr.ind=TRUE)
    mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
    mtx
}

remove_all_na <- function(betas) {
    namtx <- is.na(betas)
    remove <- rowSums(namtx) > ncol(betas) - 2
    if (sum(remove) != 0) {
        betas <- betas[-which(remove),]
        message(sprintf("removing %i rows (too many missing)",sum(remove)))
    }
    betas
}

prepare_graph_mtx <- function(betas,graphSize=10000) {
    if (nrow(betas) < graphSize) {
        graphSize <- nrow(betas)
    }
    sds <- order(-apply(betas,1,sd,na.rm=TRUE))[seq(graphSize)]
    graph_mtx <- imputeRowMean(betas[sds,])
    graph_mtx
}

build_imputation_graph <- function(graph_mtx,graphSize=10000) {
    g <- rnndescent::nnd_knn(
        graph_mtx,
        k = 50,
        metric="correlation"
    )
    rownames(g$idx) <- rownames(graph_mtx)
    rownames(g$dist) <- rownames(graph_mtx)
    g
}

prepare_imputation_search_graph <- function(graph_mtx,graph) {
    rnndescent::prepare_search_graph(
        graph_mtx,
        graph,
        metric="correlation",
        verbose = FALSE
    )
}

query_imputation_graph <- function(betas,qry,graph_mtx,searchGraph,k=2) {
    qry_mtx <- imputeRowMean(betas[qry,])
    query_nn <- rnndescent::graph_knn_query(
        query = qry_mtx,
        reference = graph_mtx,
        reference_graph = searchGraph,
        k = k,
        metric = "correlation",
        verbose = FALSE
    )
    rownames(query_nn$idx) <- rownames(qry_mtx)
    rownames(query_nn$dist) <- rownames(qry_mtx)
    query_nn
}


getDelta <- function(betas,nn_df) {
    vapply(seq(nrow(nn_df)), function(x) {
        c1 <- mean(betas[nn_df[x,1],],na.rm=TRUE)
        c2 <- mean(betas[nn_df[x,2],],na.rm=TRUE)
        c2 - c1
    },numeric(1))
}


impute <- function(mtx,nn_df,corr=0.8,delta=0.3) {
    namtx <- is.na(mtx)
    missing <- names(which(rowSums(namtx) >= 1))
    m2 <- do.call(rbind, lapply(rownames(mtx), function(x) {
        if (!x %in% missing) {
            return(mtx[x,])
        }
        nbr <- nn_df[x,]
        v <- mtx[x,]
        na <- which(is.na(v))
        if (nbr[["corr"]] >= corr & nbr[["delta"]] < delta) {
            v[na] <- mtx[nbr[["nbr"]],na]
        } else {
            v[na] <- mean(mtx[x,],na.rm=TRUE)
        }
        return(v)
    }))
    m2 <- imputeRowMean(m2)
    rownames(m2) <- rownames(mtx)
    m2
}


#' imputeMissingProbes imputes missing values based on co-methylated probes
#'
#' @param betas matrix of beta values where probes are on the
#' rows and samples on columns
#' @param graphSize # of probes to use for search graph construction
#' @param corr correlation threshold to determine imputation
#' strategy. If correlation with neighbor is below threshold, row mean
#' imputation is performed
#' @param delta mean difference threshold between imputation target
#' and neighbor. If mean difference exceeds threshold, row mean
#' imputation is performed
#' @return A betas matrix with missing values imputed
#' @examples
#' library(SummarizedExperiment)
#' se <- sesameDataGet('MM285.467.SE.tissue20Kprobes')
#' betas <- assay(se)
#' betas <- imputeMissingProbes(betas)
#' sesameDataGet_resetEnv()
#'
#' @export
imputeMissingProbes <- function(betas,graphSize=10000,corr=0.8,delta=0.3) {
    b2 <- remove_all_na(betas)
    namtx <- is.na(b2)
    graph_mtx <- prepare_graph_mtx(b2,graphSize = graphSize)
    g <- build_imputation_graph(graph_mtx = graph_mtx)
    sg <- prepare_imputation_search_graph(graph_mtx = graph_mtx,graph=g)
    qry <- names(which(rowSums(namtx) > 0))
    query_nn <- query_imputation_graph(
        betas=b2,
        qry=qry,
        graph_mtx = graph_mtx,
        searchGraph = sg
    )
    nbr_col_lgl <- rownames(query_nn$idx) == rownames(g$idx)[query_nn$idx[,1]]
    nbr_col <- ifelse(nbr_col_lgl,2,1)
    nbrs <- vapply(seq(nbr_col), function(x) {
        rownames(g$idx)[query_nn$idx[x,nbr_col[x]]]
    },character(1))
    nbr_corr <- 1 - vapply(seq(nbr_col), function(x) {
        query_nn$dist[x,nbr_col[x]]
    },numeric(1))
    nn_df <- data.frame(
        qry=qry,
        nbr=nbrs,
        corr=nbr_corr
    )
    delta_beta <- getDelta(b2,nn_df)
    nn_df[["delta"]] <- delta_beta
    rownames(nn_df) <- qry
    impute(b2,nn_df=nn_df,corr=corr,delta=delta)
}
