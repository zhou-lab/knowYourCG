
#' test whether query is significantly enriched in database
testEnrichment = function(setQ, setD, setU) {
    mtx = matrix(c(
        length(intersect(setQ, setD)),
        length(setdiff(setD, setQ)), 
        length(setdiff(setQ, setD)), 
        length(setdiff(setU, union(setD, setQ)))),
        nrow = 2, 
        dimnames = list(
            Query = c("Q_in","Q_out"),
            Database = c("D_in","D_out")))
    list(mtx = mtx, test = fisher.test(mtx))
}

#' calculate fold change given a matrix
calcFoldChange = function(mtx){
	num = mtx[1, 1] / (mtx[1, 1] + mtx[1, 2])
	den = (mtx[1, 1] + mtx[2, 1]) / sum(mtx)
	num / den
}

#' create a volcano plot given a filename, fold change (log2), pvalue (-log10), and title (optional), xlabel (optional), ylabel (optional)
plotVolcano = function(filename, fc, pvalue, title="", xlabel=NA, ylabel=NA) {

    data = na.omit(data.frame(fc=fc, 
        pvalue=pvalue))
    rownames(data) = 1:nrow(data)
    colnames(data) = c("log2fc", "pvalue")

    title = gsub('(.{1,80})(\\s|$)', '\\1\n', title)

    if (is.na(xlabel)) {
        xlabel = "log2 fold change"
    }

    if (is.na(ylabel)) {
        ylabel = "-log10 pvalue"
    }

    pdf(filename, height=50, width=50, onefile=FALSE)
    ggplot(data=data, aes(x=fc, y=pvalue, color = cut(pvalue, c(2.995732, Inf) ))) + geom_point() + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    labs(title = title, fill = "pvalue") + 
    theme(plot.title = element_text(size=80, face = "bold"), 
        axis.text=element_text(size=48), 
        axis.title=element_text(size=60), 
        legend.title = element_text(size=48),
        legend.text=element_text(size=40)) # + 
    # scale_color_manual(name = "pvalue", values = c("[2.995732, Inf)" = "red"), labels = c("> 2.995732"))
    dev.off(0)
}

#' test all databaseSet and return a list ranked by enrichment (odds-ratio)
testEnrichmentAll = function() {
	
}


