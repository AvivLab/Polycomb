# Written by Maryl Lambros on Oct 7, 2019 
# Script to analyze combined computational model results to use
# for doing similar analyses as single-cell RNA-seq datasets, like PCA


# User inputs:
folderToGetResults = "20191212-123317" #"20191107-165139" #"20191024-090520"
prcRemovedFolder = "PrcRemoved1and2"
workingDir = paste("/home/maryl/Desktop/Polycomb/julia1.1.1_MarylNewVersionWithTwoEnvs_Aug2019/dataOutput/",folderToGetResults,"/",prcRemovedFolder,"/", sep="")
setwd(workingDir)
input_file = "combinedModelData"

conditionsToPlotString = "env1ReplaceEnv2" # three possible inputs: "all", "env1ReplaceEnv2", or "env2ReplaceEnv1"
plotNormTestBool = TRUE # if want to plot count depth relationship after normalization for model sequencing-like data


# Start:
library("Seurat")
library("RColorBrewer")

data=read.delim(input_file, header=F, stringsAsFactors=F, sep="\t")
labels=data[1,]
numRows<-dim(data)[1]
geneSymbols<- paste("gene",c(1:(numRows-1)),sep="")

levelsAll<-c("AfterEvolEnv1","AfterEvolEnv2","PcgIntactEnv1","PcgIntactEnv2",
          "PcgBrokenEnv1","PcgBrokenEnv2","PcgIntactEnv1ReplaceEnv2",
          "PcgIntactEnv2ReplaceEnv1","PcgBrokenEnv1ReplaceEnv2","PcgBrokenEnv2ReplaceEnv1",
          "PcgBrokenEnv3ReplaceEnv1","PcgIntactEnv3ReplaceEnv1",
          "PcgBrokenEnv3ReplaceEnv2","PcgIntactEnv3ReplaceEnv2")

## Make white color transparent so can keep same PCA plot structure while add in conditions for presentations:
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

## Determine conditions to plot:
# All conditions: 
conditionsToPlotAll = c("PcgIntactEnv1","PcgIntactEnv2","PcgBrokenEnv2","PcgBrokenEnv1","PcgBrokenEnv2ReplaceEnv1",
                        "PcgBrokenEnv1ReplaceEnv2","PcgIntactEnv2ReplaceEnv1","PcgIntactEnv1ReplaceEnv2")
# Conditions for environment 1 replacing environment 2:
conditionsToPlotEnv1ReplaceEnv2 = c("PcgIntactEnv1","PcgIntactEnv2","PcgBrokenEnv1ReplaceEnv2","PcgIntactEnv1ReplaceEnv2")
# Conditions for env 2 replacing env 1:
conditionsToPlotEnv2ReplaceEnv1 = c("PcgIntactEnv1","PcgIntactEnv2","PcgBrokenEnv2ReplaceEnv1","PcgIntactEnv2ReplaceEnv1")
# Set which conditions to plot based on user input above:
if (conditionsToPlotString == "all") {
  conditionsToPlot = conditionsToPlotAll
  # set point style for each condition plotting:
  pts<- c(1,1,16,16,1,1,6,1)
  # set colors to use for each condition plotting:
  colorsToUse = brewer.pal(n = length(conditionsToPlot), name = "Set1")
} else if (conditionsToPlotString == "env1ReplaceEnv2") {
  conditionsToPlot = conditionsToPlotEnv1ReplaceEnv2
  # set point style for each condition plotting:
  pts<- c(1,1,1,1)
  # set colors to use for each condition plotting:
  colorsToUse = c("#E41A1C", "#377EB8", "#FF7F00", "#A65628") # to make points transparent: c("#E41A1C", "#377EB8", t_col("#FFFFFF", 100), t_col("#FFFFFF", 100)) #original colors: c("#E41A1C" (red), "#377EB8" (blue), "#FF7F00" (orange), "#A65628" (brown))
} else if (conditionsToPlotString == "env2ReplaceEnv1") {
  conditionsToPlot = conditionsToPlotEnv2ReplaceEnv1
  # set point style for each condition plotting:
  pts<- c(1,1,1,1)
  # set colors to use for each condition plotting:
  colorsToUse = c("#E41A1C", "#377EB8", "#FFFF33", "#F781BF")
}

myinds <- labels %in% conditionsToPlot
levels <- levelsAll[levelsAll %in% conditionsToPlot]
titleForPlot = "ForEnvironment1"

mylabels<-labels[myinds]
mydata<-sapply(data[2:numRows,myinds],as.numeric) # 2:numRows because model experiment type labels (represented by each column) are in 1st row
rownames(mydata)<-geneSymbols
colnames(mydata)<-mylabels
mylabels<-colnames(mydata)
f <- factor(mylabels,levels)

# pull out constants from run to get number of independent trails and population size:
constantsData <- read.delim("../constants.jl", header=T, stringsAsFactors=F)
library("stringr")
# number indepedent trails:
numTrials <- as.double(str_extract(constantsData[1,1], "[1-500]+"))
# population size:
numInds <- as.double(str_extract(constantsData[7,1], "\\d+"))


# # histograms of unnormalized data:
# dev.new()
# hist(mydata[,mylabels=="PcgBrokenEnv3ReplaceEnv2"])
# dev.new()
# hist(mydata[,mylabels=="PcgIntactEnv3ReplaceEnv1"], plot = TRUE)
# 
# normalizedMyData <- scale(log(mydata+1),center=TRUE,scale=TRUE)
# rownames(normalizedMyData)<-geneSymbols
# colnames(normalizedMyData)<-mylabels
# # Plot Histograms
# dev.new()
# hist(normalizedMyData[,mylabels=="PcgBrokenEnv3ReplaceEnv1"], plot = TRUE)
# dev.new()
# hist(normalizedMyData[,mylabels=="PcgIntactEnv3ReplaceEnv1"], plot = TRUE)
# # --> bimodal signatures are different...why not separating in PCA?
# if (plotNormTestBool) {
#   pdf(paste("CountDepthPlot",titleForPlot,"InDataFolder",folderToGetResults,".pdf",sep=""), height=5, width=7)
#   plotresults <- plotCountDepth(normalizedMyData, Conditions = mylabels, FilterExpression = 0) # test if need to do single-cell normalization or not (do values line up with 1 or curves not lining up on x axis? If not lining up, then need to use scnorm)
#   # Since gene expression increases proportionally with sequencing depth, we expect to find the
#   # estimated count-depth relationships near 1 for all genes. This is typically true for bulk RNAseq datasets. However, it does not hold in most single-cell RNA-seq datasets. In this example
#   # data, the relationship is quite variable across genes.
#   dev.off()
#   str(plotresults)
# }
# 
# pca<-prcomp(normalizedMyData,center=TRUE,scale.=TRUE)
# #mykmeans<-kmeans(pca2,2)
# #pdf(paste("PCA",titleForPlot,"InDataFolder",folderToGetResults,".pdf",sep=""), height=5, width=7)
# pts<- rep(1, length(unique(mylabels)))
# pts[grep("Broken",levels)] <- 16
# dev.new()
# plot(pca$rotation[,1],pca$rotation[,2],col=f, pch=pts[f],xlab="Principal Component 1",ylab="Principal Component 2", main = "")
# legend("bottomleft",NULL,levels(f),col=1:length(f),pch=pts)
# #dev.off()

# create Seurat object of data:
seuratObjModelData <- CreateSeuratObject(counts = mydata, min.cells = 0, min.features = 0)
seuratObjModelData <- NormalizeData(seuratObjModelData)
seuratObjModelData <- FindVariableFeatures(seuratObjModelData)
seuratObjModelData <- ScaleData(seuratObjModelData, features = geneSymbols, do.scale = FALSE, do.center = FALSE)
seuratObjModelData$CellType <- mylabels
Idents(seuratObjModelData) <- "CellType"

# dev.new()
# FeatureScatter(seuratObjModelData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.new()
VlnPlot(seuratObjModelData, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)



### Plot PCA for a single independent trial:
seuratData = seuratObjModelData@assays$RNA@scale.data
indTrialNum = 1
firstIndividualInIndependentTrial = seq((indTrialNum-1)*numInds+1,(numInds*numTrials)*length(levels),by=numInds*numTrials)
lastIndividualInIndependentTrial = seq((indTrialNum-1)*numInds+numInds,numInds*numTrials*(length(levels)),by=numInds*numTrials)
individualsInPopForGivenLevels = unlist(Map(':',firstIndividualInIndependentTrial,lastIndividualInIndependentTrial))
singlePopData = seuratData[,individualsInPopForGivenLevels]
singlePopLabels = mylabels[individualsInPopForGivenLevels]
factorSinglePop <- factor(singlePopLabels,levels)
colnames(singlePopData) <- singlePopLabels
pcaSinglePop<-prcomp(singlePopData,center=FALSE,scale=FALSE)
dev.new()
plot(pcaSinglePop$rotation[,1],pcaSinglePop$rotation[,2],col=colorsToUse[factorSinglePop], pch=pts[factorSinglePop], cex=0.9, xlab="Principal Component 1",ylab="Principal Component 2")
legend("bottomright",NULL,levels(factorSinglePop),col=colorsToUse,pch=pts)
dim(singlePopData[,singlePopLabels=="PcgIntactEnv1"])



### Plot PCA results from using seuratObjModelData but using R's prcomp function so can control colors and shapes better
seuratData = seuratObjModelData@assays$RNA@scale.data
# Set the column names back to mylabels
colnames(seuratData) <- mylabels
pca<-prcomp(seuratData,center=FALSE,scale=FALSE)
#pdf("modelResultsForEnv1ReplaceEnv2BlueRedOrangeConditionsOnlyPCA.pdf", height=5, width=7)
dev.new()
plot(pca$rotation[,1],pca$rotation[,2],col=colorsToUse[f], pch=pts[f], cex=0.7, xlab="Principal Component 1",ylab="Principal Component 2")
# add centroids of clusters:
# x <- pca$rotation[,1]
# y <- pca$rotation[,2]
# class <- mylabels
# df <- data.frame(class, x, y)
# centroids <- aggregate(cbind(x,y)~class,df,mean)
# points(centroids$x, centroids$y, col=c("yellow","green4","blue","darkred","darkblue"), pch=c(19,19,19,19,19), cex=2)
legend("topright",NULL,levels(f),col=colorsToUse,pch=rep(16,length(pts)))
#dev.off()


## Use Seurat package to plot PCA:
seuratObjModelData <- RunPCA(seuratObjModelData, features = geneSymbols, verbose = FALSE, npcs = 20)
# Plot PCA:
dev.new()
DimPlot(seuratObjModelData, reduction = "pca", cols = c("purple","green","cyan","red","blue"), group.by = "ident", pt.size = 0.3)
#title(main = paste("PCA for",titleToAdd,"when use",matchingTitle,scalingTitle,"the data for genes DE btw primary & metastatic cells"))
# Plot information about PCA:
dev.new()
VizDimLoadings(seuratObjModelData, dims = 1:2, reduction = "pca")
dev.new()
ElbowPlot(seuratObjModelData)

# Plot UMAP:
seuratObjModelData <- RunUMAP(seuratObjModelData, dims = 1:3)
dev.new()
DimPlot(seuratObjModelData, reduction = "umap", cols = c("purple","green","cyan","red","blue"), group.by = "ident", pt.size = 0.3)


## Plot PCA results for just PRC intact in env 1 vs env 2:
indsTwoDifEnvs <- labels=="PcgIntactEnv1" | labels=="PcgIntactEnv2" #| labels=="PcgIntactEnv3ReplaceEnv1" | labels=="PcgIntactEnv3ReplaceEnv2"
myLabelsTwoDifEnvs<-labels[indsTwoDifEnvs]
levels<-c("AfterEvolEnv1","AfterEvolEnv2","PcgIntactEnv1","PcgIntactEnv2",
          "PcgBrokenEnv1","PcgBrokenEnv2","PcgIntactEnv1ReplaceEnv2",
          "PcgIntactEnv2ReplaceEnv1","PcgBrokenEnv1ReplaceEnv2","PcgBrokenEnv2ReplaceEnv1",
          "PcgBrokenEnv3ReplaceEnv1","PcgIntactEnv3ReplaceEnv1",
          "PcgBrokenEnv3ReplaceEnv2","PcgIntactEnv3ReplaceEnv2")
levelIndsTwoDifEnvs<- levels=="PcgIntactEnv1" | levels=="PcgIntactEnv2" #| levels=="PcgIntactEnv3ReplaceEnv1" | levels=="PcgIntactEnv3ReplaceEnv2" 
levelsTwoDifEnvs<-levels[levelIndsTwoDifEnvs]
fTwoDifEnvs <- factor(myLabelsTwoDifEnvs,levelsTwoDifEnvs)
dataTwoDifEnvsOnly <- sapply(data[2:numRows,indsTwoDifEnvs],as.numeric) # 2:numRows because model experiment type labels (represented by each column) are in 1st row
rownames(dataTwoDifEnvsOnly)<-geneSymbols
colnames(dataTwoDifEnvsOnly)<-myLabelsTwoDifEnvs
seuratObjModelDataTwoDifEnvs <- CreateSeuratObject(counts = dataTwoDifEnvsOnly, min.cells = 0, min.features = 0)
seuratObjModelDataTwoDifEnvs <- NormalizeData(seuratObjModelDataTwoDifEnvs)
seuratObjModelDataTwoDifEnvs <- FindVariableFeatures(seuratObjModelDataTwoDifEnvs)
seuratObjModelDataTwoDifEnvs <- ScaleData(seuratObjModelDataTwoDifEnvs, features = geneSymbols, do.scale = FALSE, do.center = FALSE)
seuratObjModelDataTwoDifEnvs$CellType <- myLabelsTwoDifEnvs
Idents(seuratObjModelDataTwoDifEnvs) <- "CellType"
seuratObjModelDataTwoDifEnvs <- RunPCA(seuratObjModelDataTwoDifEnvs, features = geneSymbols, verbose = FALSE, npcs = 20)
seuratObjModelDataTwoDifEnvs <- RunUMAP(seuratObjModelDataTwoDifEnvs, dims = 1:3)
dev.new()
DimPlot(seuratObjModelDataTwoDifEnvs, reduction = "umap", cols = c("red","green"), group.by = "ident", pt.size = 0.3)
# Plot pca:
seuratDataTwoDifEnvs = seuratObjModelDataTwoDifEnvs@assays$RNA@scale.data
colnames(seuratDataTwoDifEnvs) <- myLabelsTwoDifEnvs
pts<- c(16,16,16,16)
colorsToUse = c("purple","green","cyan","red")
pcaTwoDifEnvs<-prcomp(seuratDataTwoDifEnvs,center=FALSE,scale=FALSE)
#pdf("modelResultsPCAForEnv1And2.pdf", height=5, width=7)
library("plot3D")
dev.new()
plot(pcaTwoDifEnvs$rotation[,1],pcaTwoDifEnvs$rotation[,2],col=colorsToUse[fTwoDifEnvs], pch=pts[fTwoDifEnvs], cex=0.5)
scatter3D(pcaTwoDifEnvs$rotation[,1],pcaTwoDifEnvs$rotation[,2],pcaTwoDifEnvs$rotation[,3],col=colorsToUse[fTwoDifEnvs], pch=pts[fTwoDifEnvs], cex=0.5, xlab="Principal Component 1",ylab="Principal Component 2", angle = 90)
# add centroids of clusters:
#x <- pcaTwoDifEnvs$rotation[,1]
#y <- pcaTwoDifEnvs$rotation[,2]
#class <- myLabelsTwoDifEnvs
#df <- data.frame(class, x, y)
#centroids <- aggregate(cbind(x,y)~class,df,mean)
#points(centroids$x, centroids$y, col=c("yellow","green4","blue","darkred"), pch=pts+3, cex=2)
#legend("bottomleft",NULL,levels(fTwoDifEnvs),col=colorsToUse,pch=pts)
#dev.off()



### Find PcG target genes:
polycombStateMatrix <- read.delim("polycombStateAfterEvolEnv2.csv", header= TRUE, stringsAsFactors=F, sep=",")
polycombStateAverages <- rowMeans(polycombStateMatrix)

polycombVecMatrix <- read.delim("polycombVecAfterEvolEnv2.csv", header= TRUE, stringsAsFactors=F, sep=",")
polycombVecAverages <- rowMeans(polycombVecMatrix)
  
  
### Gene by gene differences analysis for intact PcG in
# new env compared to intact PcG in env 2, and broken
# PcG in new env compared to intact PcG in env 2:
datm<-seuratData[,mylabels=='PcgBrokenEnv3ReplaceEnv2']
datp<-seuratData[,mylabels=='PcgIntactEnv3ReplaceEnv2']
datn<-seuratData[,mylabels=='AfterEvolEnv2']
meansm<-rowMeans(datm)
meansp<-rowMeans(datp)
meansn<-rowMeans(datn)
# Aim 2 H2 quantify distance results:
cat("Comp. Model: Distance btw intact in new env. & env 2:",dist(rbind(meansp,meansn)),
    "\nDistance btw broken in new env & env 2:",dist(rbind(meansm,meansn)),"\n")
# Plot each gene by gene differences results:
diff1<-meansp-meansn
diff2<-meansm-meansn
colorsGeneDiffs <- rep("black", length(diff2))
colorsGeneDiffs[polycombStateAverages <= 0.20] <- "red"
#pdf(paste("CompModelGeneByGeneDifferences.pdf",sep=""), height=5, width=7)
dev.new()
plot(diff1,diff2,main = "Gene by Gene Differences", xlab = "Distance btw intact in new env & env 2", ylab = "Distance btw broken in new env and env 2", col=colorsGeneDiffs, xlim=c(-3,3), ylim=c(-3,3))#xlim=c(-max(abs(diff1)),max(abs(diff1))), ylim=c(-max(abs(diff2)),max(abs(diff2))))
abline(h=0)
abline(v=0)
#dev.off()
c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
print(paste("Comp Model: # of genes in red =",length(colorsGeneDiffs[colorsGeneDiffs=="red"]), "out of", length(colorsGeneDiffs)))
cat(length(geneSymbols[polycombStateAverages <= 0.20]),"PRC suppressed genes:",geneSymbols[polycombStateAverages <= 0.20])
cat(paste("average log fold change for population:",mean(abs(diff2/diff1))))
cat(paste("Distance btw broken in new env and evolved env for population is",dist(rbind(meansm,meansn))))
cat(paste("Distance btw intact in new env and evolved env for population is",dist(rbind(meansp,meansn))))


### Above analysis but for only 1 given independent trial because know which are suppressed by Polycomb and which aren't better because all in same population:
numIndTrials <- 20
differenceResults <- rep(0,numIndTrials)
distanceResultsBroken <- rep(0,numIndTrials)
distanceResultsIntact <- rep(0,numIndTrials)
for (i in 1:numIndTrials) {
  indTrialIndividualsInds <- grep(paste("Trial",i,".",sep=""),colnames(polycombStateMatrix), fixed = TRUE)
  numIndividualsInPop <- length(indTrialIndividualsInds)
  datm<-seuratData[,mylabels=='PcgBrokenEnv3ReplaceEnv2'][,indTrialIndividualsInds]
  datp<-seuratData[,mylabels=='PcgIntactEnv3ReplaceEnv2'][,indTrialIndividualsInds]
  datn<-seuratData[,mylabels=='AfterEvolEnv2'][,indTrialIndividualsInds]
  meansm<-rowMeans(datm)
  meansp<-rowMeans(datp)
  meansn<-rowMeans(datn)
  # Plot each gene by gene differences results:
  diff1<-meansp-meansn
  diff2<-meansm-meansn
  colorsGeneDiffs <- rep("black", length(diff2))
  polycombStateAverages <- rowMeans(polycombStateMatrix[,indTrialIndividualsInds])
  polycombVecAverages <- rowMeans(polycombVecMatrix[,indTrialIndividualsInds])
  colorsGeneDiffs[polycombVecAverages >= 0.95] <- "cyan"
  colorsGeneDiffs[polycombStateAverages == 0.0] <- "red"
  #pdf(paste("CompModelGeneByGeneDifferences.pdf",sep=""), height=5, width=7)
  # dev.new()
  # plot(diff1,diff2,main = paste("Gene by Gene Differences for Ind. Trial", i), xlab = "Distance btw intact in new env & env 2", ylab = "Distance btw broken in new env and env 2", col=colorsGeneDiffs, xlim=c(-max(abs(diff1)),max(abs(diff1))), ylim=c(-max(abs(diff2)),max(abs(diff2))))
  # abline(h=0)
  # abline(v=0)
  # #dev.off()
  c(sqrt(sum(diff1^2)),sqrt(sum(diff2^2)),mean(abs(diff2)<abs(diff1)))
  # Aim 2 H2 quantify distance results:
  differenceResults[i] <- mean(abs(diff2/diff1))
  cat(paste("\naverage log fold change for individual",i,"is",differenceResults[i]))
  distanceResultsBroken[i] <- dist(rbind(meansm,meansn))
  distanceResultsIntact[i] <- dist(rbind(meansp,meansn))
  cat(paste("\nDistance btw broken in new env and evolved env for independent trial",i,"is",distanceResultsBroken[i]))
  cat(paste("\nDistance btw intact in new env and evolved env for independent trial",i,"is",distanceResultsIntact[i]))
  print(paste("\nComp Model: # of genes in red =",length(colorsGeneDiffs[colorsGeneDiffs=="red"]), "out of", length(colorsGeneDiffs)))
  cat(length(geneSymbols[polycombStateAverages == 0.0]),"PRC suppressed genes:",geneSymbols[polycombStateAverages == 0.0])
}

cat(paste("\nMean of distance btw broken in new env and evolved env for all independent trials =",mean(distanceResultsBroken)))
cat(paste("\nMean of distance btw intact in new env and evolved env for all independent trials =",mean(distanceResultsIntact)))
t.test(x = distanceResultsBroken, y = distanceResultsIntact, paired = TRUE)


### Look at per individual*
numIndTrials <- 20
popSize <- 500
differenceResults <- rep(0,numIndTrials)
distanceResultsBroken <- rep(0,numIndTrials*popSize)
distanceResultsIntact <- rep(0,numIndTrials*popSize)
distanceResultsBrokenInEvolvedEnv <- rep(0,numIndTrials*popSize)
for (i in 1:numIndTrials) {
  indTrialIndividualsInds <- grep(paste("Trial",i,".",sep=""),colnames(polycombStateMatrix), fixed = TRUE)
  numIndividualsInPop <- length(indTrialIndividualsInds)
  datm<-seuratData[,mylabels=='PcgBrokenEnv3ReplaceEnv2'][,indTrialIndividualsInds]
  datp<-seuratData[,mylabels=='PcgIntactEnv3ReplaceEnv2'][,indTrialIndividualsInds]
  datn<-seuratData[,mylabels=='AfterEvolEnv2'][,indTrialIndividualsInds]
  date <- seuratData[,mylabels=='PcgBrokenEnv2'][,indTrialIndividualsInds]
  popSize <- dim(datn)[2]
  for (j in 1:popSize) {
    meansm<-datm[,j]
    meansp<-datp[,j]
    meansn<-datn[,j]
    meanse <- date[,j]
    # Plot each gene by gene differences results:
    # diff1<-meansp-meansn
    # diff2<-meansm-meansn
    # colorsGeneDiffs <- rep("black", length(diff2))
    # polycombStateAverages <- rowMeans(polycombStateMatrix[,indTrialIndividualsInds])
    # polycombVecAverages <- rowMeans(polycombVecMatrix[,indTrialIndividualsInds])
    # colorsGeneDiffs[polycombVecAverages >= 0.95] <- "cyan"
    # colorsGeneDiffs[polycombStateAverages == 0.0] <- "red"
    #pdf(paste("CompModelGeneByGeneDifferences.pdf",sep=""), height=5, width=7)
    # dev.new()
    # plot(diff1,diff2,main = paste("Gene by Gene Differences for Ind. Trial", i), xlab = "Distance btw intact in new env & env 2", ylab = "Distance btw broken in new env and env 2", col=colorsGeneDiffs, xlim=c(-max(abs(diff1)),max(abs(diff1))), ylim=c(-max(abs(diff2)),max(abs(diff2))))
    # abline(h=0)
    # abline(v=0)
    # #dev.off()
    # Aim 2 H2 quantify distance results:
    # differenceResults[i] <- mean(abs(diff2/diff1))
    # cat(paste("\naverage log fold change for individual",i,"is",differenceResults[i]))
    distanceResultsBroken[j+(i-1)*popSize] <- dist(rbind(meansm,meansn))
    distanceResultsIntact[j+(i-1)*popSize] <- dist(rbind(meansp,meansn))
    distanceResultsBrokenInEvolvedEnv[j+(i-1)*popSize] <- dist(rbind(meanse,meansn))
    # cat(paste("\nDistance btw broken in new env and evolved env for independent trial",i,"individual",j,"is",distanceResultsBroken[j+(i-1)*popSize]))
    # cat(paste("\nDistance btw intact in new env and evolved env for independent trial",i,"individual",j,"is",distanceResultsIntact[j+(i-1)*popSize]))
    # print(paste("\nComp Model: # of genes in red =",length(colorsGeneDiffs[colorsGeneDiffs=="red"]), "out of", length(colorsGeneDiffs)))
    # cat(length(geneSymbols[polycombStateAverages == 0.0]),"PRC suppressed genes:",geneSymbols[polycombStateAverages == 0.0])
  }
}
cat(paste("\nMean of distance btw broken in new env and evolved env for all independent trials =",mean(distanceResultsBroken)))
cat(paste("\nMean of distance btw intact in new env and evolved env for all independent trials =",mean(distanceResultsIntact)))
t.test(x = distanceResultsBroken, y = distanceResultsIntact, paired = TRUE)
cat(paste("\nMean of distance btw broken in evolved env and evolved env for all independent trials =",mean(distanceResultsBrokenInEvolvedEnv)))
t.test(x = distanceResultsBroken, y = distanceResultsBrokenInEvolvedEnv, paired = TRUE)


### Subhypothesis 3 analysis for model data except look at pcg target genes,
# so expect that as broken pcg mechanism in new env moves farther from intact
# in original env that the genes suppressed by PcG mechanism will increase in
# expression (so positive correlation):
numIndTrials <- 10
for (i in 1:numIndTrials) {
  indTrialIndividualsInds <- grep(paste("Trial",i,".",sep=""),colnames(polycombStateMatrix), fixed = TRUE)
  numIndividualsInPop <- length(indTrialIndividualsInds)
  polycombStateAverages <- rowMeans(polycombStateMatrix[,indTrialIndividualsInds])
  numPcgTargetsGenes <- length(geneSymbols[polycombStateAverages == 0.0])
  allGenesExceptPcgMechanism <- geneSymbols[polycombStateAverages > 0.0]
  mydata2m<-seuratData[,mylabels=='PcgBrokenEnv3ReplaceEnv2'][,indTrialIndividualsInds]
  mydata2p<-seuratData[,mylabels=='PcgIntactEnv3ReplaceEnv2'][,indTrialIndividualsInds]
  mydata2n<-seuratData[,mylabels=='AfterEvolEnv2'][,indTrialIndividualsInds]
  
  len<-length(allGenesExceptPcgMechanism)#len<-sum(diffnotpcginds)
  myscores<-rep(0,dim(mydata2m)[2])
  meansn<-rowMeans(mydata2n)
  for (i in 1:dim(mydata2m)[2]) { # loop through each metastatic cell to find the distance between each cell compared to mean of normal
    metastaticCellExpression <- as.numeric(mydata2m[,i])
    myscores[i] <- dist(rbind(metastaticCellExpression,meansn))
  }
  # Overall distances:
  mymeanscores<-myscores
  
  pcgdat<-mydata2m[geneSymbols[polycombStateAverages == 0.0],]
  pcgdat<-data.frame(cbind(t(pcgdat),mymeanscores))
  colnames(pcgdat)<-c(geneSymbols[polycombStateAverages == 0.0],"score")
  
  mymodel<-lm(score~.,data=pcgdat)
  summary(mymodel)
  # Pull out PcG mechanism genes that are significant:
  pvalsOfCoefficients <- summary(mymodel)$coefficients[-1,4] # first element is the p-value for the intercept and rest are for the slopes of each PcG mechanism gene, so remove first element
  significantCoefficientInds <- pvalsOfCoefficients <= 0.1
  # PcG genes with positive coefficients, meaning
  # decreased expression positively correlated
  # with movement towards normal:
  pcgMechanismSlopes <- mymodel$coefficients[-1] # remove first element that is the intercept estimate because we just want to look at estimates of the slope for each PcG mechanism gene
  pcgMechanismSlopesSignificant <- pcgMechanismSlopes[significantCoefficientInds]# pull out only slope values that are significant
  pcgGenesPos<-pcgMechanismSlopesSignificant[pcgMechanismSlopesSignificant > 0]
  #pcgGenesPos<- mymodel$coefficients[mymodel$coefficients > 0]
  cat("# PcG genes w/ positive coefficients:",length(pcgGenesPos),
      "out of",length(pcgMechanismSlopesSignificant),"with significant coefficients")
  # Bootstrapping to find significance of results:
  controlPos<-sample(geneSymbols, numPcgTargetsGenes) # If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: controlPos<-sample(which(!pcgOnlyInds & !diffinds), numPcgMechanismOnlyGenes)
  controldat<-seuratData[,mylabels=='PcgBrokenEnv3ReplaceEnv2'][controlPos,indTrialIndividualsInds]
  controldat<-data.frame(cbind(t(controldat),mymeanscores))
  colnames(controldat)<-c(controlPos,"score") # If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: colnames(controldat)<-c(names(controlPos),"score")
  mymodel<-lm(score~.,data=controldat)
  
  nposdist<-rep(0,1000)
  for (i in 1:1000)
  {
    controlPos<-sample(geneSymbols, numPcgTargetsGenes) # If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: controlPos<-sample(which(!pcgOnlyInds & !diffinds), numPcgMechanismOnlyGenes)
    controldat<-seuratData[,mylabels=='PcgBrokenEnv3ReplaceEnv2'][controlPos,indTrialIndividualsInds]
    controldat<-data.frame(cbind(t(controldat),mymeanscores))
    colnames(controldat)<-c(controlPos,"score") # If using all genes except differentially expressed btw primary and metastatic and also pcg mechanism genes: colnames(controldat)<-c(names(controlPos),"score")
    controlmodel<-lm(score~.,data=controldat)
    npositive<-sum(controlmodel$coefficients[2:(numPcgTargetsGenes+1)][summary(controlmodel)$coefficients[-1,4] <= 0.1]>0)
    #npositive<-sum(controlmodel$coefficients[-1]>0)
    nposdist[i]<-npositive
  }
  dev.new()
  hist(nposdist, breaks=13,ylim=c(0,100))
  abline(v=length(pcgGenesPos),col="red")
  
  cat(paste("Mean and std of positive distribution:", mean(nposdist[!is.na(nposdist)]), sd(nposdist[!is.na(nposdist)])))
  library("fitdistrplus")
  fit <- fitdist(nposdist[!is.na(nposdist)], "norm")
  pValHist <- pnorm(q = length(pcgGenesPos), mean= fit$estimate[['mean']], sd= fit$estimate[['sd']])
  cat(paste("P-value of # PcG genes w/ positive coefficients:", 1-pValHist))
}




### Differential expression results
difExpBetweenTypes <- function(seuratObjData,geneSymbols,testToUse,ident1String,ident2String){
  differences <- FindMarkers(seuratObjData, ident.1 = ident1String, ident.2 = ident2String, test.use = testToUse, logfc.threshold = 0.1)
  diffinds <- sapply(geneSymbols[1:length(geneSymbols)],function(x) {any(rownames(differences)==x)})
  difExpGenes <- geneSymbols[diffinds]
  return(difExpGenes)
}
library("MAST")
difExpressedGenes <- difExpBetweenTypes(seuratObjModelData,geneSymbols,"wilcox","PcgBrokenEnv3ReplaceEnv2","PcgIntactEnv3ReplaceEnv2")
dev.new()
RidgePlot(seuratObjModelData, features = difExpressedGenes, sort = TRUE)
dev.new()
VlnPlot(seuratObjModelData, features = difExpressedGenes, sort = TRUE)


# Plot UMAP:
seuratObjModelData <- RunUMAP(seuratObjModelData, dims = 1:3)
dev.new()
DimPlot(seuratObjModelData, reduction = "umap")


# Plot TSNE
seuratObjModelData <- RunTSNE(seuratObjModelData, dims = 1:2, method = "FIt-SNE")
dev.new()
DimPlot(seuratObjModelData, reduction = "tsne")


# Differential expression analysis using Limma (no genes found to be differentially expressed 10-15-19):
library("limma")
logData=voom(mydata,plot=FALSE)
design <- model.matrix(~0+f)
colnames(design)=levels
fit=lmFit(logData,design)
logData<-data.frame(logData)
# Differential expression btw metastatic and primary cancer cells
contrastTumor=makeContrasts(ccDif=PcgBrokenEnv3ReplaceEnv1-PcgIntactEnv3ReplaceEnv1, levels=design)
fitTumor=contrasts.fit(fit, contrastTumor)
fit2=eBayes(fitTumor)
fdr=5
results <- decideTests(fit2,method="separate",lfc=1,p.value=fdr)
diffinds <- results[,1]!=0
diffPositions=which(diffinds) # get the row position in the original data set of genes with differential expression
# Data matrix containing only those genes that differentailly
# expressed btw metastatic and primary cancer cells:
differences=logData[diffPositions, ] # from the original data table, extract the rows for the differentially-expressed genes
cat("Dim. of [genes x cells] of metastatic vs primary differenital expression:", dim(differences))