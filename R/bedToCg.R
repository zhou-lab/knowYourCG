#' Convert BED CpG set to YAME .cg format
#'
#' Converts a BED file listing CpGs into a YAME-compressed .cg file using
#' the YAME reference coordinate and command-line tools.
#'
#' @param bed_file Character string. Path to BED file listing CpGs.
#' @param ref_cr Character string. Path to YAME reference coordinate file (.cr).
#' @param out_file Character string. Path to output .cg file.
#' @param sort_bed Logical. Whether to sort the BED file with bedtools before
#'   intersecting. (Default: TRUE)
#' @param verbose Logical. Whether to print the shell commands. (Default: FALSE)
#'
#' @return Invisibly returns the output file path.
#' @export
#'
#' @examples
#' # Download YAME reference coordinate file for mm10
#' ref_cr <- tempfile(fileext = ".cr")
#' download.file(
#'   "https://zenodo.org/records/18176270/files/mm10chr16_f7_cpg.cr",
#'   destfile = ref_cr, quiet = TRUE, mode = "wb"
#' )
#'
#' # Create a BED file with CpG positions
#' bed_file <- tempfile(fileext = ".bed")
#' bed_data <- c(
#'     "chr16\t3003648\t3003650",
#'     "chr16\t3003766\t3003768",
#'     "chr16\t3004036\t3004038",
#'     "chr16\t3004052\t3004054",
#'     "chr16\t3004275\t3004277",
#'     "chr16\t3004470\t3004472",
#'     "chr16\t3004545\t3004547",
#'     "chr16\t3004633\t3004635",
#'     "chr16\t3004652\t3004654",
#'     "chr16\t3004739\t3004741",
#'     "chr16\t3004841\t3004843",
#'     "chr16\t3004959\t3004961",
#'     "chr16\t3005065\t3005067",
#'     "chr16\t3005444\t3005446",
#'     "chr16\t3005522\t3005524")
#' 
#' writeLines(bed_data, bed_file)
#'
#' # Convert BED to YAME .cg format (requires bedtools and yame)
#' out_file <- tempfile(fileext = ".cg")
#' \donttest{
#' bedToCg(bed_file, ref_cr, out_file, verbose = TRUE)
#' }
#'
#' # Clean up
#' unlink(c(bed_file, ref_cr, out_file))
bedToCg <- function(bed_file, ref_cr, out_file, sort_bed = TRUE, verbose = FALSE) {
    if (!is.character(bed_file) || length(bed_file) != 1) {
        stop("'bed_file' must be a single character string.", call. = FALSE)
    }
    if (!is.character(ref_cr) || length(ref_cr) != 1) {
        stop("'ref_cr' must be a single character string.", call. = FALSE)
    }
    if (!is.character(out_file) || length(out_file) != 1) {
        stop("'out_file' must be a single character string.", call. = FALSE)
    }
    if (!file.exists(bed_file)) {
        stop("BED file not found: ", bed_file, call. = FALSE)
    }
    if (!file.exists(ref_cr)) {
        stop("Reference coordinate file not found: ", ref_cr, call. = FALSE)
    }

    bedtools_path <- Sys.which("bedtools")
    yame_path <- Sys.which("yame")
    if (bedtools_path == "") {
        stop("bedtools not found on PATH. Please install bedtools.", call. = FALSE)
    }
    if (yame_path == "") {
        stop("yame not found on PATH. Please install yame.", call. = FALSE)
    }

    bed_file <- normalizePath(bed_file, mustWork = TRUE)
    ref_cr <- normalizePath(ref_cr, mustWork = TRUE)
    out_file <- normalizePath(out_file, mustWork = FALSE)

    sorted_bed <- bed_file
    temp_bed <- NULL
    if (isTRUE(sort_bed)) {
        temp_bed <- tempfile(pattern = "kycg_bed_", fileext = ".bed")
        if (isTRUE(verbose)) {
            message(sprintf("bedtools sort -i %s > %s", bed_file, temp_bed))
        }
        status <- system2("bedtools",
            args = c("sort", "-i", bed_file),
            stdout = temp_bed,
            stderr = "")
        if (!identical(status, 0L)) {
            stop("bedtools sort failed. Check your BED file.", call. = FALSE)
        }
        sorted_bed <- temp_bed
    }

    cmd_pack <- sprintf(
        "yame unpack %s | bedtools intersect -a - -b %s -c -sorted | cut -f4 | yame pack -fb - > %s",
        shQuote(ref_cr),
        shQuote(sorted_bed),
        shQuote(out_file)
    )
    if (isTRUE(verbose)) message(cmd_pack)
    status <- system2("sh", args = c("-c", cmd_pack), stdout = FALSE, stderr = "")
    if (!identical(status, 0L)) {
        stop("YAME/bedtools pipeline failed. Check input files and tools.", call. = FALSE)
    }

    if (!is.null(temp_bed) && file.exists(temp_bed)) {
        unlink(temp_bed)
    }

    invisible(out_file)
}
