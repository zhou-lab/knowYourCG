#' plotVolcano creates a volcano plot of -log2(p.value) and log(estimate)
#' given data with fields estimate and p.value.
#'
#' @param data DataFrame where each field is a database name with two fields
#' for the estimate and p.value.
#' @param title String representing the title label. Optional. (Default: NA)
#' @param xlabel String representing the x-axis label. Optional. (Default: 
#' NA)
#' @param ylabel String representing the y-axis label. Optional. (Default: 
#' NA)
#'
#' @return A DataFrame with the estimate/statistic, p-value, and name of test 
#' for the given results.
#'
#' @export
plotVolcano = function(data, title=NA, xlabel=NA, ylabel=NA) {
    if (is.na(title)) {
        title = "Volcano plot"
    } 
    title = gsub('(.{1,80})(\\s|$)', '\\1\n', title)

    data = data.frame(estimate=log(data$estimate, base=10), p.value=-log(data$p.value, base=2))

    if (is.na(xlabel)) {
       xlabel = "log2 fold change"
    }

    if (is.na(ylabel)) {
       ylabel = "-log10 pvalue"
    }

    # pdf(filename, height=50, width=50, onefile=FALSE)
    if (any(data$p.value > -log(0.05))) {
        g = ggplot(data=data, aes(x=estimate, y=p.value, color = cut(p.value, c(-log(0.05), Inf) ))) 
    } else {
        g = ggplot(data=data, aes(x=estimate, y=p.value)) 
    }
    g = g + geom_point() + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    labs(title = title, legend= "HELLO", fill = "pvalue") + 
    theme(plot.title = element_text(size=16, face = "bold"), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=12), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
    scale_colour_discrete(name = "Significance (alpha = 0.05)", 
        labels=c("Significant", "Not Significant"))
        #labs(color='Legend Title') # + 
    # scale_color_manual(name = "pvalue", values = c("[2.995732, Inf)" = "red"), labels = c("> 2.995732"))
    # dev.off(0)
    g
}

