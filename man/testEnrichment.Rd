% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testEnrichment.R
\name{testEnrichment}
\alias{testEnrichment}
\title{testEnrichment tests for the enrichment of set of probes (query set) in
a number of features (database sets).}
\usage{
testEnrichment(
  query,
  databases = NULL,
  universe = NULL,
  alternative = "greater",
  platform = NULL,
  silent = FALSE
)
}
\arguments{
\item{query}{Vector of probes of interest (e.g., significant probes)}

\item{databases}{List of vectors corresponding to the database sets of
interest with associated meta data as an attribute to each element. Optional.
(Default: NA)}

\item{universe}{Vector of probes in the universe set containing all of
the probes to be considered in the test. If it is not provided, it will be
inferred from the provided platform. (Default: NA).}

\item{alternative}{"two.sided", "greater", or "less"}

\item{platform}{String corresponding to the type of platform to use. Either
MM285, EPIC, HM450, or HM27. If it is not provided, it will be inferred
from the query set probeIDs (Default: NA).}

\item{silent}{output message? (Default: FALSE)}
}
\value{
One list containing features corresponding the test estimate,
p-value, and type of test.
}
\description{
testEnrichment tests for the enrichment of set of probes (query set) in
a number of features (database sets).
}
\examples{

library(SummarizedExperiment)
df <- rowData(sesameDataGet('MM285.tissueSignature'))
query <- df$Probe_ID[df$branch == "B_cell"]
res <- testEnrichment(query, "chromHMM", platform="MM285")
sesameDataGet_resetEnv()

}