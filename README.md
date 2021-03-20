# KnowYourCpG - DNAme tool

NEED TO UPDATE

To install from github {Not currently implemented},
```R
install.packages("BiocManager")
library(BiocManager)
install("ethanmoyer/knowyourcpg")
```

Just clone to install for now.

## Install KnowYourCpG

Development version can be installed from github
```{r, eval=FALSE}
BiocManager::install('ethanmoyer/knowyourcpg')
BiocManager::install('ethanmoyer/knowyourcpgData')
```

```{r}
target_set = c("cg12345","cg12346")
testEnrichment(target_set, reference = "PMDs", background = "EPIC")
scanAllBiology(target_set, background = "EPIC")
## ranked list of enrichment ordered by odds ratio
##--------------------
## PMDs 5-fold enrichment, pval=0.00012
## repeat 4-fold enrichment, pval=0.004

```

## Background set
- Infinium EPIC probes
- HM450 probes
- All CpGs in the human genome

## Reference set

### technical
- Probes with measurement artifact
- Infinium-I vs Infinium-II
- cg vs ch vs rs
- cg content
- meQTLs (methylation linked to genetic mutations)

### biology
- Ancestry-associated CpGs
- specific genes
- Polycomb target (H3K27me3 marked)
- PMDs
- solo-WCGWs
- CpG island
- promoter CpGs
- CTCF CpGs
- repetitive elements
- genes on specific pathway






