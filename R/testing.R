rm(list=ls())

source('~/projects/tbmate/scripts/tbmate.R')
source('~/projects/knowYourCpG/R/analyze.R')

mm285index <- '~/Dropbox/Ongoing_knowYourCpG/TBK_INDICES/MM285.idx.gz'

InfiniumType <- c("/Users/Fred/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285//20210409_Infinium_Type.tbk")
infiniumDB <- tbk_data(tbk_fnames = InfiniumType, idx_fname = mm285index)

# universal probe set
probeIDs <- names(infiniumDB)

# fake vectors of significant probe names
# should return no associations
randomList <- sample(names(infiniumDB), 600)

# should be enriched for type I probes
infiniumList <- c(
  sample(names(infiniumDB)[infiniumDB == 'I'], 3000), 
  sample(names(infiniumDB)[infiniumDB == 'II'], 2)
)

cpgdensity <- c("/Users/Fred/Dropbox/Ongoing_knowYourCpG/DATABASE_SETS/MM285/20210416_cpg_density.tbk")
cpgdensityDB <- tbk_data(tbk_fnames = cpgdensity, idx_fname = mm285index)

# should be enriched for high cpg density
densityList <- sample(names(sort(cpgdensityDB, decreasing = TRUE))[1:1000], 300)

# generate test results
randomTest <- testEnrichmentAll(probeIDs = probeIDs, sigProbes = randomList)
print(randomTest$categorical$`20210409_Infinium_Type`$I$Fisher.Test)
print(randomTest$continuous)

infiniumTest <- testEnrichmentAll(probeIDs = probeIDs, sigProbes = infiniumList)
print(infiniumTest$categorical$`20210409_Infinium_Type`$I$Fisher.Test)
print(infiniumTest$categorical$`20210409_Infinium_Type`$II$Fisher.Test)
print(infiniumTest$continuous)

densityTest <- testEnrichmentAll(probeIDs = probeIDs, sigProbes = densityList)
print(densityTest$categorical$`20210409_Infinium_Type`$I$Fisher.Test)
print(densityTest$continuous)
