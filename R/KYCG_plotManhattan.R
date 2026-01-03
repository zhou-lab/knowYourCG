#' KYCG_plotManhattan makes a manhattan plot to summarize EWAS results
#'
#' @param vals named vector of values (P,Q etc), vector name is Probe ID.
#' @param platform String corresponding to the type of platform to use for
#' retrieving GRanges coordinates 
#' of probes. Either MM285, EPIC, HM450, or HM27. If it is not provided, it
#' will be inferred from the query set probeIDs (Default: NA).
#' @param genome hg38, mm10, ..., will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @param title title for plot
#' @param label_min Threshold above which data points will be labelled with
#' Probe ID
#' @param col color
#' @param ylabel y-axis label
#' @param rasterize if true use ggrastr to rasterize non-significant data.
#' @param rasterize_thres, the threshold of rasterize
#' @return a ggplot object
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges end
#' @importFrom tibble as_tibble
#' @import ggrepel
#' @import sesameData
#' @examples
#'
#' library(sesameData)
#' 
#' ## Create example with simulated -log10(p-values)
#' ## Mix of non-significant (low values) and significant (high values)
#' probes <- names(sesameData_getManifestGRanges("HM450"))
#' set.seed(123)
#' vals <- setNames(
#'     c(runif(990, 0, 3),      # Non-significant probes
#'       runif(10, 5, 25)),     # Significant probes
#'     sample(probes, 1000)
#' )
#'
#' KYCG_plotManhattan(vals,
#'     platform = "HM450",
#'     title = "Example Manhattan Plot",
#'     ylabel = "-log10(P-value)",
#'     label_min = 20)
#'
#' @export
KYCG_plotManhattan <- function(
    vals, platform = NULL, genome = NULL, title = NULL,
    rasterize=FALSE, rasterize_thres = 3,
    label_min = 100, col = c("wheat1", "sienna3"), ylabel="Value") {

    stopifnot(is(vals, "numeric"))
    if (is.null(platform)) { platform <- queryCheckPlatform(
        platform, query=vals, silent = FALSE) }
    genome <- sesameData_check_genome(genome, platform)
    gr <- sesameData_getManifestGRanges(platform, genome=genome)
    seqLength <- sesameData_getGenomeInfo(genome)$seqLength
    
    v <- vals[names(gr)]
    gr <- gr[!is.na(v)]
    SummarizedExperiment::mcols(gr)$val <- v[!is.na(v)]
    cumLength <- setNames(c(0,cumsum(
        as.numeric(seqLength))[-length(seqLength)]), names(seqLength))
    midLength <- cumLength + seqLength/2
    SummarizedExperiment::mcols(gr)$pos <- cumLength[
        as.character(seqnames(gr))] + GenomicRanges::end(gr)
    SummarizedExperiment::mcols(gr)$Probe_ID <- names(gr)
    
    df <- as_tibble(gr)
    df$seqnames <- factor(df$seqnames, levels=names(seqLength))
    requireNamespace("ggrepel")
    p <- ggplot(df, aes(x=.data[["pos"]], y=.data[["val"]],
        color=.data[["seqnames"]]))
    if (rasterize) {
        requireNamespace("ggrastr")
        p <- p + ggrastr::rasterise(geom_point(
            alpha=0.5, size=1,
            data = df[df$val < rasterize_thres,]), dpi=600)
    } else {
        p <- p + geom_point(
            alpha=0.5, size=1,
            data = df[df$val < rasterize_thres,])
    }
    p <- p + geom_point(
        alpha=0.8, size=1,
        data = df[df$val >= rasterize_thres,])
    p <- p + ggrepel::geom_text_repel(data=df[df$val > label_min,],
        aes(label = .data[["Probe_ID"]])) +
        scale_color_manual(values = rep(col, length(seqLength))) +
        scale_x_continuous(labels = names(midLength), breaks= midLength) +
        scale_y_continuous(expand = c(0, 0)) +  
        theme_bw() +
        theme(
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) + labs(title=title) + xlab("Chromosome") + ylab(ylabel)
    p
}

