## Offline tests for knowledgebase loading, annotation, per-set statistics,
## probe proximity and the aggregate plots. Synthetic data only.

## getDBs() names each set by its "dbname" attribute, not by its list name,
## so synthetic sets carry both, as a loaded knowledgebase does.
as_dbs2 <- function(dbs, group = "synthetic") {
    setNames(lapply(names(dbs), function(nm) {
        d <- dbs[[nm]]
        attr(d, "group") <- group
        attr(d, "dbname") <- nm
        d
    }), names(dbs))
}

write_kb <- function(path, kb) {
    tbl <- data.frame(Probe_ID = unlist(kb, use.names = FALSE),
                      Knowledgebase = rep(names(kb), lengths(kb)))
    readr::write_tsv(tbl, path)
}

test_that("loadDBs reads a directory or a vector of files", {
    d <- file.path(tempdir(), "kycg_loaddbs")
    dir.create(d, showWarnings = FALSE)
    write_kb(file.path(d, "G1.gz"), list(a = c("cg1", "cg2"), b = "cg3"))
    write_kb(file.path(d, "G2.gz"), list(c = c("cg4", "cg5", "cg6")))

    dbs <- loadDBs(d)
    expect_setequal(names(dbs), c("a", "b", "c"))
    expect_equal(attr(dbs$a, "group"), "G1")        # .gz dropped
    expect_equal(attr(dbs$c, "dbname"), "c")
    expect_equal(as.vector(dbs$c), c("cg4", "cg5", "cg6"))

    one <- loadDBs(file.path(d, "G2.gz"))
    expect_equal(names(one), "c")

    expect_setequal(listDBGroups(path = d), c("G1.gz", "G2.gz"))
    unlink(d, recursive = TRUE)
})

test_that("listDBGroups lists the built-in catalogue, filtered", {
    gps <- listDBGroups()
    expect_true(all(c("Title", "Description") %in% colnames(gps)))
    expect_gt(nrow(gps), 10)
    mm <- listDBGroups("MM285")
    expect_true(all(grepl("MM285", mm$Title)))
    expect_lt(nrow(mm), nrow(gps))
})

test_that("annoProbes labels probes by set membership", {
    dbs <- list(A = c("cg1", "cg2"), B = c("cg2", "cg3"))
    probes <- c("cg1", "cg2", "cg3", "cg4")
    a <- annoProbes(probes, dbs, platform = "EPIC", silent = TRUE)
    expect_equal(unname(a), c("A", "A,B", "B", NA))
    expect_equal(names(a), probes)

    a2 <- annoProbes(probes, dbs, platform = "EPIC", sep = ";",
                     silent = TRUE)
    expect_equal(unname(a2[2]), "A;B")

    m <- annoProbes(probes, dbs, platform = "EPIC", indicator = TRUE,
                    silent = TRUE)
    expect_true(is.matrix(m))
    expect_equal(dim(m), c(4, 2))
    expect_equal(colnames(m), c("A", "B"))
    expect_equal(unname(m["cg2", ]), c(TRUE, TRUE))

    ## unnamed sets are named db1, db2, ...
    u <- annoProbes(probes, unname(dbs), platform = "EPIC", indicator = TRUE,
                    silent = TRUE)
    expect_equal(colnames(u), c("db1", "db2"))

    ## no databases: all NA
    n <- annoProbes(probes, NULL, platform = "EPIC", silent = TRUE)
    expect_true(all(is.na(n)))
    expect_equal(names(n), probes)
})

test_that("dbStats summarizes a beta matrix per set", {
    set.seed(1)
    ids <- sprintf("cg%04d", 1:100)
    betas <- matrix(runif(300), 100, 3, dimnames = list(ids, c("s1", "s2", "s3")))
    betas[1:5, "s3"] <- NA
    dbs <- as_dbs2(list(first = ids[1:10], last = ids[91:100],
                        absent = c("cgX", "cgY")))

    st <- dbStats(betas, dbs)
    expect_equal(dim(st), c(3, 3))
    expect_equal(rownames(st), c("s1", "s2", "s3"))
    expect_equal(colnames(st), c("first", "last", "absent"))
    expect_equal(st["s1", "first"], mean(betas[1:10, "s1"]))
    expect_equal(st["s3", "first"], mean(betas[6:10, "s3"]))  # NA dropped
    expect_true(all(is.na(st[, "absent"])))

    ## n_min above the non-missing count masks the statistic
    st2 <- dbStats(betas, dbs, n_min = 8)
    expect_true(is.na(st2["s3", "first"]))
    expect_false(is.na(st2["s1", "first"]))

    ## another function, a plain vector, and the long form
    med <- dbStats(betas, dbs, fun = median)
    expect_equal(med["s2", "last"], median(betas[91:100, "s2"]))
    v <- dbStats(betas[, "s1"], dbs)
    expect_equal(rownames(v), "sample")
    lg <- dbStats(betas, dbs, long = TRUE)
    expect_true(all(c("query", "db", "value") %in% colnames(lg)))
    expect_equal(nrow(lg), 9)

    ## unnamed sets fall back to the dbname attribute
    st3 <- dbStats(betas, unname(dbs))
    expect_equal(colnames(st3), c("first", "last", "absent"))
})

test_that("testProbeProximity finds clustered probes on a given GRanges", {
    ## 200 probes on one chromosome, 10 kb apart, except a run of five
    ## packed 100 bp apart
    pos <- seq(1, by = 10000, length.out = 200)
    pos[101:105] <- pos[100] + (1:5) * 100
    ids <- sprintf("cg%04d", 1:200)
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, width = 2))
    names(gr) <- ids

    set.seed(1)
    res <- testProbeProximity(ids[c(100:105, 1, 50, 150)], gr = gr,
                              iterations = 20)
    expect_named(res, c("Stats", "Clusters"))
    expect_equal(res$Stats$num_query, 9)
    expect_equal(res$Stats$hits_query, 5)
    expect_lt(res$Stats$p.val, 0.01)
    expect_true(all(c("seqnames", "start", "end", "distance") %in%
                    colnames(res$Clusters)))
    expect_equal(nrow(res$Clusters), 6)

    ## nothing within the bin: p = 1, no clusters
    far <- testProbeProximity(ids[c(1, 50, 150)], gr = gr, iterations = 5)
    expect_equal(far$Stats$hits_query, 0)
    expect_equal(far$Stats$p.val, 1)
    expect_true(is.na(far$Clusters))

    ## a null with no hits at all is forced to a small lambda
    set.seed(2)
    expect_warning(lone <- testProbeProximity(ids[c(100, 101)], gr = gr,
                                              iterations = 5),
                   "Forcing lambda")
    expect_equal(lone$Stats$lambda, 0.001)
})

test_that("aggregate plots accept a list of enrichment results", {
    set.seed(3)
    mk <- function() data.frame(dbname = paste0("st", 1:6),
                                estimate = rnorm(6))
    rl <- list(s1 = mk(), s2 = mk(), s3 = mk())
    expect_s3_class(KYCG_plotPointRange(rl), "ggplot")

    meta <- function() data.frame(dbname = as.character(1:5),
                                  db = 1:5,
                                  label = paste0("bin", 1:5),
                                  estimate = rnorm(5))
    expect_s3_class(KYCG_plotMetaEnrichment(list(a = meta(), b = meta())),
                    "ggplot")
    expect_s3_class(KYCG_plotMetaEnrichment(meta()), "ggplot")
    expect_error(KYCG_plotMetaEnrichment(list(mk())))
})

test_that("KYCG_plotEnrichAll draws a multi-group result", {
    set.seed(4)
    n <- 40
    df <- data.frame(
        dbname = paste0("db", seq_len(n), ";short", seq_len(n)),
        group = rep(c("KYCG.MM285.chromHMM.20210210",
                      "KYCG.MM285.TFBSconsensus.20220116",
                      "KYCG.MM285.gene.00000000",
                      "KYCG.MM285.seqContext.20210630"), each = 10),
        estimate = abs(rnorm(n)) + 0.5,
        FDR = 10^-runif(n, 3, 40))
    df$gene_name <- paste0("Gene", seq_len(n))
    pe <- function(...) expect_no_warning(KYCG_plotEnrichAll(...))
    expect_s3_class(pe(df), "ggplot")
    expect_s3_class(pe(df, short_label = FALSE, n_label = 3), "ggplot")

    plain <- df[, c("dbname", "group", "estimate", "FDR")]
    plain$group <- rep(c("G1", "G2"), each = 20)
    expect_s3_class(pe(plain), "ggplot")
})

test_that("bedToCg checks its arguments before running anything", {
    bed <- tempfile(fileext = ".bed")
    writeLines("chr1\t0\t10", bed)
    expect_error(bedToCg(1, bed, "o.cg"), "'bed_file' must be")
    expect_error(bedToCg(bed, 1, "o.cg"), "'ref_cr' must be")
    expect_error(bedToCg(bed, bed, 1), "'out_file' must be")
    expect_error(bedToCg("no_such.bed", bed, "o.cg"), "BED file not found")
    expect_error(bedToCg(bed, "no_such.cr", "o.cg"),
                 "Reference coordinate file not found")
    unlink(bed)
})

test_that("bedToCg turns a BED into a format 0 .cg", {
    skip_on_os("windows")
    skip_if(!nzchar(Sys.which("bedtools")), "bedtools not on PATH")
    skip_if(!nzchar(Sys.which("yame")), "yame not on PATH")
    td <- tempfile("kycg_bed2cg_"); dir.create(td)
    ## ten CpGs on chr1, one per 100 bp, as a reference coordinate file
    cr_txt <- file.path(td, "ref.bed")
    writeLines(sprintf("chr1\t%d\t%d\tcg%02d", (0:9) * 100,
                       (0:9) * 100 + 2, 1:10), cr_txt)
    cr <- file.path(td, "ref.cr")
    system2("yame", c("pack", "-f", "r", cr_txt, cr))
    bed <- file.path(td, "regions.bed")
    writeLines(c("chr1\t250\t450", "chr1\t0\t50"), bed)
    out <- file.path(td, "out.cg")
    expect_message(bedToCg(bed, cr, out, verbose = TRUE), "bedtools")
    expect_true(file.exists(out))
    txt <- system2("yame", c("unpack", out), stdout = TRUE)
    expect_equal(as.integer(txt), c(1, 0, 0, 1, 1, 0, 0, 0, 0, 0))
    unlink(td, recursive = TRUE)
})
