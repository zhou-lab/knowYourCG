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

## bit-wise masking quick calculation

## Universe set
- Infinium EPIC probes - 850k
- HM450 probes
- All CpGs in the human genome - WGBS 28M

## Database set

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
- TFBS
- chromatin accessibility, tumor biology
- PMDs
- solo-WCGWs
- CpG island
- promoter CpGs
- CTCF CpGs
- repetitive elements
- genes on specific pathway
- sequence conservation
- imprinting
- X/Y
- tissue-specific methylation signature
- age-associated methylation

### Query sets

## Notes

- Look for manifest files here 

https://zwdzwd.github.io/InfiniumAnnotation


