#' databaseSetGet retrieves database sets from a meta data sheet by querying the array, group, and reference columns. The data is returned as a list where the names correspond to chosen database sets.
#'
#' @param meta List for the meta data corresponding to the information about each database set.
#' @param array Vector of array platforms for the query.
#' database sets. Optional. (Default: NULL)
#' @param group String corresponding to the type of platform to use. Either 
#' MM285, EPIC, HM450, or HM27.
#' @param reference Logical value indicating whether to display intermediate 
#' text output about the type of test. Optional. (Default: FALSE)
#'
#' @return One list containing features corresponding the test estimate, p-value, and type of test.
#'
#' @export
databaseSetGet = function(meta, array=c("EPIC", "HM450", "MM285"), group=c("PMDsoloWCGW", "CpGIsland"), reference=NA) {

	meta_ = meta

	if (!is.na(array)) {
		meta = meta[meta$Array %in% array, ]
	}

	if (!is.na(group)) {
		meta = meta[meta$Group %in% group, ]
	}

	if (!is.na(reference)){
		meta = meta[meta$Reference %in% reference]
	}

	data = flattenlist(lapply(meta$Location, function(x) {readRDS(url(x))}))

	return(data)
}

flattenlist <- function(x){  
	morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
	out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
	if(sum(morelists)){ 
    	Recall(out)
 	} else{
    	return(out)
	}
}