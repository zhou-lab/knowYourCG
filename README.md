<p align="center">
<img src="https://github.com/user-attachments/assets/ee2328b9-5516-4b28-9fce-9e456534b3b1" width="220">
</p>
<p align="center">
<a href="https://github.com/zhou-lab/knowYourCG/commits/devel"><img src="https://img.shields.io/github/last-commit/zhou-lab/knowYourCG.svg?style=flat-square" alt="Last Commit"></a>
<a href="https://github.com/zhou-lab/knowYourCG/actions/workflows/build-on-windows.yaml"><img src="https://github.com/zhou-lab/knowYourCG/actions/workflows/build-on-windows.yaml/badge.svg" alt="CI Status"></a>
</p>

<!-- Release: [![Bioconductor Build Status (release)](https://bioconductor.org/shields/build/release/bioc/knowYourCG.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/knowYourCG) -->
<!-- Devel: [![Bioconductor Build Status (devel)](https://bioconductor.org/shields/build/devel/bioc/knowYourCG.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/knowYourCG) -->


# KnowYourCG (Supervised Learning of Sparse DNA Methylome Data)

## Install knowYourCG

Release version can be installed from Bioconductor
```r
BiocManager::install("knowYourCG")
```

Development version can be installed from github.
```r
BiocManager::install('zhou-lab/knowYourCG')
```

## Main Functions

### Enrichment Testing
- **`testEnrichment`** - Test for enrichment of query in knowledgebase sets
- **`testEnrichment2`** - Test enrichment from YAME-compressed CG sets
- **`testEnrichmentSEA`** - GSEA-like test for association of categorical
variable against continuous variable
- **`testProbeProximity`** - Test if query probes share closer genomic
proximity than random
- **`aggregateTestEnrichments`** - Aggregate test enrichment results

### Database Management
- **`getDBs`** - Get databases by full or partial names
- **`listDBGroups`** - List database group names
- **`loadDBs`** - Load knowledgebase databases from TSV files
- **`buildGeneDBs`** - Build gene-probe association database
- **`dbStats`** - Aggregate methylation over database set features
- **`annoProbes`** - Annotate Probe IDs using KYCG databases

### Visualization
- **`KYCG_plotBar`** - Bar plot of most enriched CG groups
- **`KYCG_plotDot`** - Dot plot of most enriched CG groups
- **`KYCG_plotEnrichAll`** - Plot enrichment test results
- **`KYCG_plotLollipop`** - Lollipop plot of log(estimate)
- **`KYCG_plotManhattan`** - Manhattan plot for EWAS results
- **`KYCG_plotMeta`** - Plot meta gene or other meta genomic features
- **`KYCG_plotMetaEnrichment`** - Plot meta gene enrichment
- **`KYCG_plotPointRange`** - Point range plot for enrichment results
- **`KYCG_plotSetEnrichment`** - Plot set enrichment
- **`KYCG_plotVolcano`** - Volcano plot of -log2(p.value) vs log(estimate)
- **`KYCG_plotWaterfall`** - Waterfall plot of log(estimate)

### Data Access
- **`kycgDataCache`** - Cache KnowYourCG data
- **`kycgDataGet`** - Get KnowYourCG data

### Utilities
- **`linkProbesToProximalGenes`** - Find genes in genomic proximity to
Infinium probes
- **`bedToCg`** - Convert BED CpG set to YAME .cg format

## See also

https://github.com/zhou-lab/YAME



