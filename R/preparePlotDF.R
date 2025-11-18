#' @importFrom dplyr slice_min
#' @importFrom dplyr ungroup
#' @importFrom stringr str_split
#' @importFrom utils head
#' @importFrom magrittr %>%
preparePlotDF <- function(
    df, n, order_by, short_label = FALSE, label_by = "dbname") {
    ## suppress R CMD CHECK no visible binding warning
    db1 <- FDR <- NULL
    
    stopifnot("estimate" %in% colnames(df) && "FDR" %in% colnames(df))
    df1 <- df[df$nD >0,]
    df1$FDR[df1$FDR==0] <- .Machine$double.xmin # cap FDR

    if (!short_label) {
        if ("group" %in% colnames(df1)) {
            gp <- sprintf("%s~", vapply(str_split(
                df1$group, "\\."), function(x) {
                    if(length(x)>3) {x[3]} else {x[1]}}, character(1)))
        } else if ("MFile" %in% colnames(df1)) {
            gp <- sprintf("%s~", sub(".cm","",df1$MFile))
        } else {
            gp <- ""
        }
    } else {
        gp <- ""
    }
    ## TODO: make everything use dbname
    if (label_by %in% colnames(df1)) {
        df1$db1 <- paste0(gp, df1[[label_by]])
    } else if ("feat" %in% colnames(df1)) { # genome-wide data
        df1$db1 <- paste0(gp, df1$feat)
    } else if ("Mask" %in% colnames(df1)) {
        df1$db1 <- paste0(gp, df1$Mask)
    }
    if (length(unique(df1$db1)) != nrow(df1)) {
        df1 <- df1 %>% group_by(db1) %>% slice_min(
            order_by=FDR, n=1, with_ties=FALSE) %>% ungroup()
    }

    ord <- df1[[order_by]]
    if (order_by == "estimate") { ord <- -ord; }
    df1 <- df1[order(ord, -df1$estimate),]
    df1 <- head(df1, n=n)

    df1$db1 <- factor(df1$db1, levels=rev(df1$db1))
    df1
}
