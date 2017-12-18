#' Filter expression matrix
#'
#' @param DataMat expression matrix (rows are cells and cols are genes)
#' @param OutThr.Gene_Expression mad on the gene expression level
#' @param OutThr.Gene_Count mad on the gene count
#' @param OutThr.Gene_Space mad on the distance from the centroid in the Gene expression space
#'
#' @return
#' @export
#'
#' @examples
FilterData <- function(DataMat,
                       OutThr.Gene_Expression = Inf,
                       OutThr.Gene_Count = Inf,
                       OutThr.Gene_Space = Inf,
                       LogSpace = FALSE,
                       PlotDebug = FALSE) {

  if(PlotDebug){
    hist(apply(DataMat > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
         freq = TRUE, ylab = "Number of cells")

    hist(apply(DataMat, 1, sum), main = "Reads per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of cells")

    hist(apply(DataMat>0, 2, sum), main = "Transcripts per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of genes identified")
  }

  # Filtering Gene count
  if(!is.infinite(OutThr.Gene_Count)){
    OutCount <- scater::isOutlier(rowSums(DataMat>0), nmads = OutThr.Gene_Count)
  } else {
    OutCount <- rep(FALSE, nrow(DataMat))
  }

  # Filtering Gene expression
  if(!is.infinite(OutThr.Gene_Expression)){
    OutExpr <- scater::isOutlier(rowSums(DataMat), nmads = OutThr.Gene_Expression)
  } else {
    OutExpr <- rep(FALSE, nrow(DataMat))
  }

  SpaceFil <- OutExpr | OutCount

  if(LogSpace){
    DataMat <- log10(DataMat + 1)
  }

  if(!is.infinite(OutThr.Gene_Space)){
    Centroid <- colMeans(DataMat[!OutExpr & !OutCount,])
    dim(Centroid) <- c(1, length(Centroid))

    DistFromCent <- as.vector(RFastDistance::fastPdist(Centroid, DataMat[!OutExpr & !OutCount,]))
    SpaceFil[!OutExpr & !OutCount] <- scater::isOutlier(DistFromCent, nmads = OutThr.Gene_Count)

  }

  return(SpaceFil)

}






#' Estimate proliferative cells
#'
#' @param DataMat expression matrix (rows are cells and cols are genes)
#' @param EstProlif estimation algorithm
#' @param QuaThr selection parameter
#'
#' @return
#' @export
#'
#' @examples
EstimateProliferativeCell <- function(DataMat, EstProlif = "MeanPerc", QuaThr = .5) {

  print("Estimating most likely proliferative cells")
  RankedData <- apply(DataMat, 1, rank)

  if(EstProlif == "Quantile"){
    print("Using quantile separation")
    NonG0Cell <-  rownames(DataMat)[apply(DataMat, 1, median) > quantile(DataMat, QuaThr)]
  }

  if(EstProlif == "PercQuant"){
    print("Using quantile ordering")
    NonG0Cell <-  rownames(DataMat)[order(apply(DataMat, 1, quantile, QuaThr), decreasing = TRUE)[1:round(nrow(DataMat)/10)]]
  }

  if(EstProlif == "KmeansPerc"){
    print("Using kmeans on quantile data")
    Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)

    KM <- kmeans(x = Cellvect, centers = range(Cellvect))
    if(PlotDebug){
      boxplot(apply(DataMat/rowSums(DataMat), 1, median) ~ KM$cluster)
    }

    if(KM$centers[1] < KM$centers[2]){
      NonG0Cell <-  rownames(DataMat)[KM$cluster == 2]
    } else {
      NonG0Cell <-  rownames(DataMat)[KM$cluster == 1]
    }

  }

  if(EstProlif == "MeanPerc"){
    print("Using mean of quantiles")
    Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
    NonG0Cell <- rownames(DataMat)[Cellvect > mean(Cellvect)]
    boxplot(Cellvect)
    abline(h=mean(Cellvect))
  }

  if(EstProlif == "AdaptiveMeanPerc"){
    print("Using adaptive mean of quantiles")
    Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
    AdjFact = 0
    while((sum(Cellvect != 0) < MinProlCells) & (QuaThr + AdjFact < 1)){
      AdjFact <- AdjFact + .01
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr + AdjFact)
    }
    boxplot(Cellvect, main=paste("Q =", QuaThr + AdjFact))
    abline(h=mean(Cellvect))
    NonG0Cell <- rownames(DataMat)[Cellvect > mean(Cellvect)]
  }

  print(paste(length(NonG0Cell), "strongly proliferative cells inferred"))

  G0Cell <-  setdiff(rownames(DataMat), NonG0Cell)

  # if(!is.null(NonG0Cell)){
  #   TB <- rbind(table(Categories[rownames(DataMat) %in% NonG0Cell]), table(Categories[rownames(DataMat) %in% G0Cell]))
  #   barplot(TB/rowSums(TB), ylab = "Percentage of cells", xlab = "Category", beside = TRUE,
  #           legend.text = c("Strongly Proliferative", "Not strongly proliferative"))
  #   barplot(TB, ylab = "Number of cells", xlab = "Category", beside = TRUE,
  #           legend.text = c("Strongly Proliferative", "Not strongly proliferative"))
  #
  # }

  return(list(G0Cell = G0Cell, NonG0Cell = NonG0Cell))

}
