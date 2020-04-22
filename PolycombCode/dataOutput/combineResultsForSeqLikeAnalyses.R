# Written by Maryl Lambros on Oct 4, 2019 
# Script to combine computational model results to use for doing similar analyses as single-cell RNA-seq datasets, like PCA

# Set folder to get results from & then set working directory to be this folder:
folderToGetResults = "20191216-105537"
prcRemovedFolder = "PrcRemoved1and2"
workingDir = paste("/home/maryl/Desktop/Polycomb/julia1.1.1_MarylNewVersionWithTwoEnvs_Aug2019/dataOutput/",folderToGetResults,"/",prcRemovedFolder,"/", sep="")
setwd(workingDir)

resultsLab <-c("AfterEvolEnv1","AfterEvolEnv2","PcgIntactEnv1","PcgIntactEnv2",
               "PcgBrokenEnv1","PcgBrokenEnv2","PcgIntactEnv1ReplaceEnv2",
               "PcgIntactEnv2ReplaceEnv1","PcgBrokenEnv1ReplaceEnv2","PcgBrokenEnv2ReplaceEnv1",
               "PcgBrokenEnv3ReplaceEnv1","PcgIntactEnv3ReplaceEnv1",
               "PcgBrokenEnv3ReplaceEnv2","PcgIntactEnv3ReplaceEnv2")
resultsTitles <- c("geneExpAfterEvolEnv1","geneExpAfterEvolEnv2","geneExpResultsPcgIntactEnv1","geneExpResultsPcgIntactEnv2",
                   "geneExpResultsPcgBrokenEnv1","geneExpResultsPcgBrokenEnv2","geneExpResultsPcgIntactEnv1ReplaceEnv2",
                   "geneExpResultsPcgIntactEnv2ReplaceEnv1","geneExpResultsPcgBrokenEnv1ReplaceEnv2","geneExpResultsPcgBrokenEnv2ReplaceEnv1",
                   "geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv1","geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv1",
                   "geneExpResultsWhenPcgBrokenAndIntroduceNewEnvInPlaceOfEnv2","geneExpResultsWhenPcgIntactAndIntroduceNewEnvInPlaceOfEnv2")
resultsFiles <- list.files(path = ".")
# remove constants.jl file from combining data sets below:
resultsFiles <- resultsFiles[!resultsFiles=="constants.jl"]
for (i in 1:length(resultsTitles)) {
  fileOfDatasetToAddToDf <- grep(paste("^",resultsTitles[i],".csv",sep=""),resultsFiles)
  # make sure only 1 file with a given result type is present; if more than one, then throw an error
  if (length(fileOfDatasetToAddToDf) == 1) {
    dataset <- read.delim(resultsFiles[fileOfDatasetToAddToDf],header = FALSE, stringsAsFactors=F, sep=",")
  } else if (length(fileOfDatasetToAddToDf) < 1) {
    print(paste("Error: no file with this name", resultsTitles[i]))
    break
  } else if (length(fileOfDatasetToAddToDf) > 1) {
    print(paste("Error: more than one file with this name", resultsTitles[i]))
    break
  }
  matrixNums <- sapply(dataset[2:dim(dataset)[1],1:dim(dataset)[2]], as.numeric)
  #matrixNums <- t(matrixNums)
  # Add result type label:
  resultLabels <- rep(resultsLab[i],dim(matrixNums)[2])
  matNames <- rbind(resultLabels,matrixNums)
  # Combine the data:
  if (i == 1) {
    allData <- matNames
  } else {
    allData <- cbind(allData,matNames)
  }
}
allData <- as.data.frame(allData)
rownames(allData) <- 1:dim(allData)[1]
colnames(allData) <- paste0("V", 1:dim(allData)[2])
gz80 = gzfile("combinedModelData", "w")
write.table(allData, gz80, sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
close(gz80)

