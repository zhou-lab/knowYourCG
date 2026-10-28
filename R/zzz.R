.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "knowYourCG: CpG-level enrichment for methylation arrays and ",
        "sequencing data.\n",
        "If it helps your work, please cite Goldberg et al., Sci Adv 2025: ",
        "citation(\"knowYourCG\")")
}
