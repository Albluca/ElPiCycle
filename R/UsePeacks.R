#' Title
#'
#' @param DataStruct
#' @param ProcStruct
#' @param MinNodes
#' @param FiltMax
#' @param Thr
#' @param QuantSel
#' @param SinglePeack
#'
#' @return
#' @export
#'
#' @examples
GetGenesWithPeaks <- function(DataStruct,
                              ProcStruct,
                              Mode = "SmoothOnCircleNodes",
                              MinNodes = 2,
                              Thr = .9,
                              QuantSel = .75,
                              Span =.3,
                              MaxPsCV = Inf,
                              nCores = 1) {


  # DataStruct = Data.Sasa.CV
  # ProcStruct = Proc.Exp.Sasa.CV
  # FiltMax = .1
  # Thr = .75
  # MinNodes = 1
  # AllGenes = TRUE
  # MedianCVThr = 1
  # QuantSel = .5


  if(Mode == "SmoothOnCircleNodes"){

    ExpMat <- DataStruct$ExpMat[, rownames(DataStruct$ExpMat) %in% names(ProcStruct$CellsPT)]
    PathLen <- ProcStruct$NodesPT[length(ProcStruct$NodesPT)]

    TopPT <- sort(ProcStruct$CellsPT, decreasing = FALSE)[1:ceiling(length(ProcStruct$CellsPT)*Span)]
    LowPT <- sort(ProcStruct$CellsPT, decreasing = TRUE)[1:ceiling(length(ProcStruct$CellsPT)*Span)]
    
    NewSpan <- Span/(1+2*Span)
    
    ExtPT <- c(
      TopPT + PathLen,
      ProcStruct$CellsPT,
      LowPT - PathLen
    )

    # Selected <- (ExtPT >= - PathLen*Span) &
    #   (ExtPT <= PathLen + PathLen*Span)

    # ExtPT <- ExtPT[Selected]
    # SelCellID <- rep(1:nrow(ExpMat), 3)[Selected]

    # Sorted <- sort(ExtPT, index.return=TRUE)
    Sorted <- sort(ExtPT)

    # SortedSelCellID <- SelCellID[Sorted$ix]

    NodesPT <- ProcStruct$NodesPT
    
    if(MaxPsCV < Inf){
      FitFun <- function(x){
        LOE <- loess(x ~ Sorted, span = NewSpan)
        predict(LOE, data.frame(Sorted = NodesPT), se = TRUE)
      }
    } else {
      FitFun <- function(x){
        LOE <- loess(x ~ Sorted, span = NewSpan)
        predict(LOE, data.frame(Sorted = NodesPT), se = FALSE)
      }
    }
    
    GeneVar <- apply(ExpMat, 2, var)
    GeneSel <- GeneVar >= quantile(GeneVar, QuantSel)

    if(nCores <= 1){

      print(paste("Computing loess smoothers on", sum(GeneSel), "genes and", length(Sorted), "pseudotime points on a single processor. This may take a while ..."))
      
      tictoc::tic()
      AllFit <- apply(ExpMat[names(Sorted),GeneSel], 2, FitFun)
      tictoc::toc()

    } else {

      no_cores <- parallel::detectCores()

      if(nCores > no_cores){
        nCores <- no_cores
        print(paste("Too many cores selected!", nCores, "will be used"))
      }

      if(nCores == no_cores){
        print("Using all the cores available. This will likely render the system unresponsive untill the operation has concluded ...")
      }

      print(paste("Computing loess smoothers on", sum(GeneSel), "genes and", length(Sorted), "pseudotime points using", nCores, "processors. This may take a while ..."))

      tictoc::tic()
      cl <- parallel::makeCluster(nCores)

      parallel::clusterExport(cl=cl, varlist=c("Sorted", "NewSpan", "NodesPT"),
                              envir = environment())

      AllFit <- parallel::parApply(cl, ExpMat[names(Sorted),GeneSel], 2, FitFun)

      parallel::stopCluster(cl)
      tictoc::toc()

    }

    if(MaxPsCV < Inf){
      PredMat <- sapply(AllFit, function(x){x$fit})
      PseudoCVMat <- sapply(AllFit, function(x){x$se.fit/x$fit})
      BinMat.CV <- (PseudoCVMat <= MaxPsCV)
    } else {
      PredMat <- AllFit
    }
    
    

    PredMat <- apply(PredMat, 2, function(x){
       x <- (x - min(x))
       x <- x/max(x)
    })


    BinMat.UP <- (PredMat >= Thr)
    BinMat.DOWN <- (PredMat <= 1-Thr)

    if(MaxPsCV < Inf){
      BinMat.UP <- BinMat.UP & BinMat.CV
      BinMat.DOWN <- BinMat.DOWN & BinMat.CV
    }

    BinMat.UP <- BinMat.UP[,colSums(BinMat.UP) > MinNodes]
    BinMat.DOWN <- BinMat.DOWN[,colSums(BinMat.DOWN) > MinNodes]

    Genes.UP <- apply(BinMat.UP, 2, which)
    Genes.DOWN <- apply(BinMat.DOWN, 2, which)

    DetMat.UP <- sapply(Genes.UP, function(x){
      sapply(ProcStruct$StageOnNodes, function(y){
        any(x %in% y)
      })
    })

    DetMat.DOWN <- sapply(Genes.DOWN, function(x){
      sapply(ProcStruct$StageOnNodes, function(y){
        any(x %in% y)
      })
    })

    UpList <- apply(DetMat.UP, 1, function(x){names(which(x))})
    DownList <- apply(DetMat.DOWN, 1, function(x){names(which(x))})

    for(PhName in union(names(UpList), names(DownList))){
      Both <- intersect(UpList[[PhName]], DownList[[PhName]])
      UpList[[PhName]] <- setdiff(UpList[[PhName]], Both)
      DownList[[PhName]] <- setdiff(DownList[[PhName]], Both)
    }

    return(list(UP = UpList, DOWN = DownList))

  }





#
#
#   if(AllGenes){
#
#     GeneExprMat.DF <- data.frame(DataStruct$ExpMat[names(ProcStruct$CellsPT), ])
#
#     DF <- PartitionData(X = DataStruct$Analysis$FinalExpMat,
#                         NodePositions = DataStruct$Analysis$FinalStruct$NodePositions,
#                         SquaredX = rowSums(DataStruct$Analysis$FinalExpMat^2))
#
#     GeneExprMat.DF.Split <- split(GeneExprMat.DF, DF$Partition)
#     GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
#     GeneExprMat.DF.Split.Sd <- lapply(GeneExprMat.DF.Split, function(x){apply(x, 2, sd)})
#
#     GeneExprMat.DF.Split.Mean.Bind <- sapply(GeneExprMat.DF.Split.Mean, cbind)
#     GeneExprMat.DF.Split.Sd.Bind <- sapply(GeneExprMat.DF.Split.Sd, cbind)
#
#     rownames(GeneExprMat.DF.Split.Mean.Bind) <- colnames(GeneExprMat.DF)
#     colnames(GeneExprMat.DF.Split.Mean.Bind) <- names(GeneExprMat.DF.Split.Mean)
#
#     rownames(GeneExprMat.DF.Split.Sd.Bind) <- colnames(GeneExprMat.DF)
#     colnames(GeneExprMat.DF.Split.Sd.Bind) <- names(GeneExprMat.DF.Split.Sd)
#
#     Reord <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
#     Reord <- Reord[Reord %in% colnames(GeneExprMat.DF.Split.Mean.Bind)]
#
#     GeneExprMat.DF.Split.Mean.Bind <- GeneExprMat.DF.Split.Mean.Bind[,Reord]
#     GeneExprMat.DF.Split.Sd.Bind <- GeneExprMat.DF.Split.Sd.Bind[,Reord]
#
#     MedianCV <- apply(GeneExprMat.DF.Split.Sd.Bind/GeneExprMat.DF.Split.Mean.Bind, 1, median, na.rm=TRUE)
#
#     # NormExp <- ProcStruct$NodesExp
#
#     NormExp <- GeneExprMat.DF.Split.Mean.Bind
#     # dim(NormExp)
#
#   } else {
#
#     NormExp <- ProcStruct$NodesExp
#     colnames(NormExp) <- ProcStruct$ExtPath
#
#     MedianCV <- rep(0, nrow(ProcStruct$NodesExp))
#
#   }
#
#
#   MedianCV <- MedianCV[rowSums(NormExp>0) > MinNodes]
#   NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
#
#   # dim(NormExp)
#
#   SDVect <- apply(NormExp, 1, sd)
#
#   SDVect[is.na(SDVect)] <- 0
#   MedianCV[is.na(MedianCV)] <- 0
#
#   NormExp <- NormExp[SDVect > quantile(SDVect, QuantSel) & MedianCV < MedianCVThr,]
#
#   NormExp[NormExp <= FiltMax] <- 0
#   NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
#
#   NormExp <- (NormExp - apply(NormExp, 1, min))/(apply(NormExp, 1, max) - apply(NormExp, 1, min))
#
#   pheatmap::pheatmap(NormExp, cluster_cols = FALSE)
#
#   FiltStagesOnNodes <- list()
#
#   for(stg in StagesNames){
#     Selected <- ProcStruct$StageOnNodes[grep(stg, names(ProcStruct$StageOnNodes))]
#     FiltStagesOnNodes[[stg]] <- unlist(Selected, use.names = FALSE)
#     names(FiltStagesOnNodes[[stg]]) <- unlist(lapply(Selected, names), use.names = FALSE)
#   }
#
#
#   StageGenes <- lapply(FiltStagesOnNodes, function(x){
#
#     SelX <- intersect(names(x), colnames(NormExp))
#     UnSelX <- setdiff(colnames(NormExp), SelX)
#
#     if(length(SelX)>0){
#       if(length(SelX)==1){
#         NormExp[, UnSelX] < Thr & NormExp[, SelX] > Thr
#       } else {
#         apply(NormExp[, UnSelX] < Thr, 1, all) & apply(NormExp[, SelX] > Thr, 1, any)
#       }
#     } else {
#       NA
#     }
#
#   })
#
#   StageGenes.Names.UP <- lapply(StageGenes, function(x){
#     names(which(x))
#   })
#
#
#   tt <- lapply(as.list(1:length(StageGenes.Names.UP)), function(i) {
#
#     if(length(StageGenes.Names.UP[[i]])>1){
#       pheatmap::pheatmap(NormExp[StageGenes.Names.UP[[i]],], cluster_cols = FALSE, main = names(StageGenes.Names.UP)[i])
#     }
#
#     if(length(StageGenes.Names.UP[[i]])==1){
#       pheatmap::pheatmap(NormExp[StageGenes.Names.UP[[i]],], cluster_cols = FALSE, cluster_rows = FALSE, main = names(StageGenes.Names.UP)[i])
#     }
#
#   })
#
#   barplot(unlist(lapply(StageGenes.Names.UP, length)), las = 2, horiz = FALSE, main = "Up genes")
#
#
#
#
#
#
#
#
#
#
#
#
#
#   StageGenes <- lapply(FiltStagesOnNodes, function(x){
#
#     SelX <- intersect(names(x), colnames(NormExp))
#     UnSelX <- setdiff(colnames(NormExp), SelX)
#
#     if(length(SelX)>0){
#       if(length(SelX)==1){
#         NormExp[, UnSelX] > 1 - Thr & NormExp[, SelX] < 1 - Thr
#       } else {
#         apply(NormExp[, UnSelX] > 1- Thr, 1, all) & apply(NormExp[, SelX] < 1 - Thr, 1, any)
#       }
#     } else {
#       NA
#     }
#
#   })
#
#   StageGenes.Names.DOWN <- lapply(StageGenes, function(x){
#     names(which(x))
#   })
#
#
#   tt <- lapply(as.list(1:length(StageGenes.Names.DOWN)), function(i) {
#
#     if(length(StageGenes.Names.DOWN[[i]])>1){
#       pheatmap::pheatmap(NormExp[StageGenes.Names.DOWN[[i]],], cluster_cols = FALSE, main = names(StageGenes.Names.DOWN)[i])
#     }
#
#     if(length(StageGenes.Names.DOWN[[i]])==1){
#       pheatmap::pheatmap(NormExp[StageGenes.Names.DOWN[[i]],], cluster_cols = FALSE, cluster_rows = FALSE, main = names(StageGenes.Names.DOWN)[i])
#     }
#
#   })
#
#   barplot(unlist(lapply(StageGenes.Names.DOWN, length)), las = 2, horiz = FALSE, main = "DOWN genes")
#
#
#
#
#
#
#
#
#
#   return(list(UP = StageGenes.Names.UP, DOWN = StageGenes.Names.DOWN))

}















#' Title
#'
#' @param PeakList
#' @param StageVect
#' @param Mode
#'
#' @return
#' @export
#'
#' @examples
ExtractStages <- function(PeakList,
                          StageVect = list(G1 = "G1", S = "S", G2M = c("G2M", "G2/M")),
                          Mode = "Any") {

  if(Mode == "Any"){

    Recoded <- lapply(PeakList, function(PhaseList) {
      lapply(StageVect, function(PhaseName){

        AllIdxs <- sapply(PhaseName, function(LocName){
          grep(pattern = LocName, x = names(PhaseList))
        })

        unlist(PhaseList[unique(unlist(AllIdxs))], use.names = FALSE)
      })
    })

    return(Recoded)

  }

  stop("Other extraction modes are not implemented yet")

}































# DataStruct = Data.Buet
# ProcStruct = Proc.Exp.Buet
# FiltMax = 0
# Thr = .7
# QuantSel = .5
# SinglePeack = TRUE
# Mode = "UP"
# MinNodes = 2
# StageInfo <- list(G1 = G1.All, S = S.All, G2M = G2M.All)
# StageInfo <- list(G1 = G1.Sel, S = S.Sel, G2M = G2M.Sel)
# StageInfo <- list(G0 = G0.All, G1 = G1.All, S = S.All, G2M = G2.All)

#' Title
#'
#' @param DataStruct
#' @param ProcStruct
#' @param StageInfo
#' @param MinNodes
#' @param FiltMax
#' @param Thr
#' @param QuantSel
#' @param Mode
#'
#' @return
#' @export
#'
#' @examples
StageWithPeaks <- function(DataStruct,
                           ProcStruct,
                           StageInfo.UP,
                           StageInfo.DOWN = list(),
                           ComputeG0 = TRUE,
                           CCGenes = NULL,
                           MinNodes = 2,
                           QuantSel = .75,
                           G0Level = .75,
                           CCLevel = .95,
                           # SinglePeack = FALSE,
                           Title = '',
                           StagingMode = "SmoothOnCircleNodes",
                           Mode = 1,
                           Span = .2,
                           nCores = 1,
                           # Thr = .9,
                           # MaxPsCV = 2,
                           FullStaging = FALSE
                              ) {

  if(is.null(CCGenes)){
    CCGenes <- unlist(StageInfo.UP, use.names = FALSE)
  }

  GenesToUse <- c(
    unlist(StageInfo.UP, use.names = FALSE),
    unlist(StageInfo.UP, use.names = FALSE),
    CCGenes
  )

  GenesToUse <- unique(GenesToUse)

  if(length(StageInfo.DOWN)==0){
    StageInfo.DOWN <- StageInfo.UP
    for(i in 1:length(StageInfo.DOWN)){
      StageInfo.DOWN[[i]] <- numeric(0)
    }
  }

  if(StagingMode == "SmoothOnCircleNodes"){

    ExpMat <- DataStruct$ExpMat[rownames(DataStruct$ExpMat) %in% names(ProcStruct$CellsPT), ]
    PathLen <- ProcStruct$NodesPT[length(ProcStruct$NodesPT)]

    GeneVar <- apply(ExpMat, 2, var)
    GeneSel <- GeneVar >= quantile(GeneVar, QuantSel)
    GeneSel <- GeneSel & (colnames(ExpMat) %in% GenesToUse)

    TopPT <- sort(ProcStruct$CellsPT, decreasing = FALSE)[1:ceiling(length(ProcStruct$CellsPT)*Span)]
    LowPT <- sort(ProcStruct$CellsPT, decreasing = TRUE)[1:ceiling(length(ProcStruct$CellsPT)*Span)]
    
    NewSpan <- Span/(1+2*Span)
    
    ExtPT <- c(
      TopPT + PathLen,
      ProcStruct$CellsPT,
      LowPT - PathLen
    )
    
    # Selected <- (ExtPT >= - PathLen*Span) &
    #   (ExtPT <= PathLen + PathLen*Span)

    # ExtPT <- ExtPT[Selected]
    # SelCellID <- rep(1:nrow(ExpMat), 3)[Selected]

    SelCellID <- rep(1:nrow(ExpMat), 3)

    Sorted <- sort(ExtPT)

    # SortedSelCellID <- SelCellID[Sorted$ix]

    NodesPT <- ProcStruct$NodesPT

    # UpdatedSpan <- length(ProcStruct$CellsPT)*Span/sum(Selected)

    FitFun <- function(x){
      LOE <- loess(x ~ Sorted, span = NewSpan)
      predict(LOE, data.frame(Sorted = NodesPT), se = FALSE)
    }
    
    # if(MaxPsCV < Inf){
    #   FitFun <- function(x){
    #     LOE <- loess(x ~ Sorted, span = NewSpan)
    #     predict(LOE, data.frame(Sorted = NodesPT), se = TRUE)
    #   }
    # } else {
    #   
    # }
    
    
    
    if(nCores <= 1){
      
      print(paste("Computing loess smoothers on", sum(GeneSel), "genes and", length(Sorted), "pseudotime points on a single processor. This may take a while ..."))
      
      op <- pboptions(type = "timer")
      PredMat <- pbapply::pbapply(ExpMat[names(Sorted),GeneSel], 2, FitFun)
      pboptions(op)
     
    } else {
      
      no_cores <- parallel::detectCores()
      
      if(nCores > no_cores){
        nCores <- no_cores
        print(paste("Too many cores selected!", nCores, "will be used"))
      }
      
      if(nCores == no_cores){
        print("Using all the cores available. This will likely render the system unresponsive untill the operation has concluded ...")
      }
      
      print(paste("Computing loess smoothers on", sum(GeneSel), "genes and", length(Sorted), "pseudotime points using", nCores, "processors. This may take a while ..."))
      
      tictoc::tic()
      cl <- parallel::makeCluster(nCores)
      
      parallel::clusterExport(cl=cl, varlist=c("Sorted", "NewSpan", "NodesPT"),
                              envir = environment())
      parallel::clusterSetRNGStream(cl, iseed = 0L)
      
      AllFit <- pbapply::pbapply(ExpMat[names(Sorted),GeneSel], 2, FitFun, cl = cl)
      # AllFit <- parallel::parApply(cl, ExpMat[names(Sorted),GeneSel], 2, FitFun)
      
      parallel::stopCluster(cl)
      tictoc::toc()
      
    }
    

    # PredMat <- sapply(AllFit, function(x){x$})
    # if(MaxPsCV < Inf){
    #   PseudoCVMat <- sapply(AllFit, function(x){x$se.fit/x$fit})
    # }

    PredMat <- apply(PredMat, 2, function(x){
      x <- (x - min(x))
      x <- x/max(x)
    })

    PredMat <- PredMat[,!apply(is.nan(PredMat), 2, any)]
    
    # BinMat.UP <- PredMat >= Thr
    # BinMat.DOWN <- PredMat <= 1-Thr
    #
    # BinMat.CV <- PseudoCVMat <= MaxPsCV
    #
    # BinMat.UP <- BinMat.UP & BinMat.CV
    # BinMat.DOWN <- BinMat.DOWN & BinMat.CV
    #
    # BinMat.UP <- BinMat.UP[,colSums(BinMat.UP) > MinNodes]
    # BinMat.DOWN <- BinMat.DOWN[,colSums(BinMat.DOWN) > MinNodes]
    #
    # Genes.UP <- apply(BinMat.UP, 2, which)
    # Genes.DOWN <- apply(BinMat.DOWN, 2, which)

    UpMat <- sapply(StageInfo.UP, function(Genes){
      if(sum(colnames(PredMat) %in% Genes)>1){
        rowSums(PredMat[, colnames(PredMat) %in% Genes])
      } else {
        rep(0, nrow(PredMat))
      }
    })


    # UpCount <- sapply(StageInfo.UP, function(Genes){
    #   sum(colnames(BinMat.UP) %in% Genes)
    # })

    DownMat <- sapply(StageInfo.DOWN, function(Genes){
      if(sum(colnames(PredMat) %in% Genes)>1){
        rowSums(1 - PredMat[, colnames(PredMat) %in% Genes])
      } else {
        rep(0, nrow(PredMat))
      }
    })

    # DownMat <- sapply(StageInfo.DOWN, function(Genes){
    #   if(sum(colnames(BinMat.DOWN) %in% Genes) > 1){
    #     rowSums(BinMat.DOWN[, colnames(BinMat.DOWN) %in% Genes])
    #   } else {
    #     rep(0, nrow(BinMat.DOWN))
    #   }
    #
    # })

    # apply(NormUp, 2, which.max)
    # apply(NormDown, 2, which.max)

    # DownCount <- sapply(StageInfo.DOWN, function(Genes){
    #   sum(colnames(BinMat.DOWN) %in% Genes)
    # })

    # NormTot[is.na(NormTot)] <- 0

    # UpMat[is.na(UpMat)] <- 0
    # DownMat[is.na(DownMat)] <- 0

    # pheatmap::pheatmap(UpMat,
    #                    cluster_rows = FALSE,
    #                    cluster_cols = FALSE)
    #
    # pheatmap::pheatmap(DownMat,
    #                    cluster_rows = FALSE,
    #                    cluster_cols = FALSE)
    #
    # pheatmap::pheatmap(t(t(UpMat+DownMat)/(DownCount + UpCount)),
    #                    cluster_rows = FALSE,
    #                    cluster_cols = FALSE)





    # Compute matrices and normalize the results


    NormDown <- t(DownMat) - apply(DownMat, 2, min)
    NormDown <- NormDown/apply(NormDown, 1, quantile, CCLevel)
    NormDown[is.na(NormDown)] <- 0
    NormDown[NormDown>1] <- 1

    NormUp <- t(UpMat) - apply(UpMat, 2, min)
    NormUp <- NormUp/apply(NormUp, 1, quantile, CCLevel)
    NormUp[is.na(NormUp)] <- 0
    NormUp[NormUp>1] <- 1

    NormTot <- t(DownMat+UpMat) - apply(DownMat+UpMat, 2, min)
    NormTot <- NormTot/apply(NormTot, 1, quantile, CCLevel)
    NormTot[is.na(NormTot)] <- 0
    NormTot[NormTot>1] <- 1


    if(ComputeG0){

      # Computing G0 (Staging matrices need to be augmented)

      MeanExp <- apply(ExpMat[, colnames(ExpMat) %in% CCGenes], 1, mean)

      SmoothMGE <- FitFun(x = MeanExp[names(Sorted)])
      SmoothMGE <- SmoothMGE-min(SmoothMGE)
      SmoothMGE <- SmoothMGE - quantile(SmoothMGE, G0Level)
      SmoothMGE[SmoothMGE > 0] <- 0
      SmoothMGE[SmoothMGE < 0] <- 1.2


      NormTot <- rbind(
        SmoothMGE,
        NormTot
      )

      NormDown <- rbind(
        SmoothMGE,
        NormDown
      )

      NormUp <- rbind(
        SmoothMGE,
        NormUp
      )

      rownames(NormTot)[1] <- "G0"
      rownames(NormUp)[1] <- "G0"
      rownames(NormDown)[1] <- "G0"

      if(any(NormDown>0)){
        pheatmap::pheatmap(NormDown, main = paste(Title, "Down"),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE)
      }

      if(any(NormUp>0)){
        pheatmap::pheatmap(NormUp, main = paste(Title, "Up"),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE)
      }

      pheatmap::pheatmap(NormTot, main = paste(Title, "Up+Down"),
                         cluster_rows = FALSE,
                         cluster_cols = FALSE)

      plot(apply(NormUp, 2, which.max), pch = 1, ylim = c(1, 4),
           yaxt = "n", ylab = "", xaxt = "n", xlab="Node ID",
           main = paste(Title, "Naive Stage assignment"))
      points(apply(NormDown, 2, which.max), pch = 2)
      points(apply(NormTot, 2, which.max), pch = 3)
      axis(side = 2, at = 1:4, labels = c("G0", "G1", "S", "G2/M"), las = 2)
      axis(side = 1, at = 1:ncol(NormDown), labels = ProcStruct$ExtPath, las = 2)
      legend(x = 1, y = 3.5, legend = c("Up", "Down", "Comb"), pch = c(1,2,3))


      print("Getting stage association")

      plot(0, 0,
           ylim = c(1, 4), xlim = c(1, ncol(NormTot)-1),
           type="n",
           yaxt = "n", ylab = "", xaxt = "n", xlab="Node ID",
           main = paste(Title, "Optimized stage assignment"))

      BestStagingTot <- AssociteNodes(PerMat = NormTot[,-ncol(NormTot)], Mode = Mode)
      for(i in 1:nrow(BestStagingTot)){
        points(BestStagingTot[i,], pch = 1)
      }

      BestStagingUp <- NULL
      BestStagingDown <- NULL

      if(FullStaging){

        if( sum(rowSums(NormUp)>0) > 1){
          BestStagingUp <- AssociteNodes(PerMat = NormUp[,-ncol(NormUp)], Mode = Mode)
          for(i in 1:nrow(BestStagingUp)){
            points(BestStagingUp[i,], pch = 2)
          }
        }

        if( sum(rowSums(NormDown)>0) > 1){
          BestStagingDown <- AssociteNodes(PerMat = NormDown[,-ncol(NormDown)], Mode = Mode)
          for(i in 1:nrow(BestStagingDown)){
            points(BestStagingDown[i,], pch = 3)
          }
        }

      }


      axis(side = 2, at = 1:4, labels = c("G0", "G1", "S", "G2/M"), las = 2)
      axis(side = 1, at = 1:(ncol(NormDown)-1),
           labels = ProcStruct$ExtPath[-length(ProcStruct$ExtPath)], las = 2)


      InfStages <- list()


      for(i in 1:nrow(BestStagingTot)){
        InfStages[[i]] <- rownames(NormTot)[BestStagingTot[i, ]]
        names(InfStages[[i]]) <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
        # InfStages[[i]] <- (InfStages[[i]])[-length(InfStages[[i]])]
      }


    } else {

      if(any(NormDown>0)){
        pheatmap::pheatmap(NormDown, main = paste(Title, "Down"),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE)
      }

      if(any(NormUp>0)){
        pheatmap::pheatmap(NormUp, main = paste(Title, "Up"),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE)
      }

      pheatmap::pheatmap(NormTot, main = paste(Title, "Up+Down"),
                         cluster_rows = FALSE,
                         cluster_cols = FALSE)

      par(mfcol = c(2,1))

      plot(apply(NormUp, 2, which.max), pch = 1, ylim = c(1, 3),
           yaxt = "n", ylab = "", xaxt = "n", xlab="Node ID",
           main = paste(Title, "Naive Stage assignment"))
      points(apply(NormDown, 2, which.max), pch = 2)
      points(apply(NormTot, 2, which.max), pch = 3)
      axis(side = 2, at = 1:3, labels = c("G1", "S", "G2/M"), las = 2)
      axis(side = 1, at = 1:ncol(NormDown), labels = ProcStruct$ExtPath, las = 2)
      legend(x = 1, y = 2.5, legend = c("Up", "Down", "Comb"), pch = c(1,2,3))

      print("Getting stage association")

      plot(0, 0,
           ylim = c(1, 3), xlim = c(1, ncol(NormTot)-1),
           type="n",
           yaxt = "n", ylab = "", xaxt = "n", xlab="Node ID",
           main = paste(Title, "Optimized stage assignment"))

      BestStagingTot <- AssociteNodes(PerMat = NormTot[,-ncol(NormTot)], Mode = Mode)
      for(i in 1:nrow(BestStagingTot)){
        points(BestStagingTot[i,], pch = 1)
      }


      BestStagingUp <- NULL
      BestStagingDown <- NULL

      if(FullStaging){

        if( sum(rowSums(NormUp)>0) > 1){
          BestStagingUp <- AssociteNodes(PerMat = NormUp[,-ncol(NormUp)], Mode = Mode)
          for(i in 1:nrow(BestStagingUp)){
            points(BestStagingUp[i,], pch = 2)
          }
        }

        if( sum(rowSums(NormDown)>0) > 1){
          BestStagingDown <- AssociteNodes(PerMat = NormDown[,-ncol(NormDown)], Mode = Mode)
          for(i in 1:nrow(BestStagingDown)){
            points(BestStagingDown[i,], pch = 3)
          }
        }

      }

      axis(side = 2, at = 1:3, labels = c("G1", "S", "G2/M"), las = 2)
      axis(side = 1, at = 1:(ncol(NormDown)-1),
                             labels = ProcStruct$ExtPath[-length(ProcStruct$ExtPath)], las = 2)

      InfStages <- list()

      for(i in 1:nrow(BestStagingTot)){
        InfStages[[i]] <- rownames(NormTot)[BestStagingTot[i, ]]
        names(InfStages[[i]]) <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
        # InfStages[[i]] <- (InfStages[[i]])[-length(InfStages[[i]])]
      }

    }



    InfStages_Ext <- InfStages
    for(i in 1:length(InfStages)){
      for(j in 1:length(InfStages[[i]])){

        PrevID <- j-1
        if(PrevID == 0){
          PrevID = length(InfStages[[i]])
        }

        NextID <- j+1
        if(NextID > length(InfStages[[i]])){
          NextID = 1
        }

        if(InfStages[[i]][j] == InfStages[[i]][PrevID] &
           InfStages[[i]][j] == InfStages[[i]][NextID]){

          InfStages_Ext[[i]][j] <- InfStages[[i]][j]

        } else {

          SelStages <- InfStages[[i]][c(PrevID, j, NextID)]
          SelStages <- SelStages[!duplicated(SelStages)]

          # print(SelStages)

          InfStages_Ext[[i]][j] <- paste(SelStages, collapse = "-")

        }

      }
    }


    # apply(NormTot, 2, which.max)
    # apply(NormTot, 2, which.max)





    return(
      list(
        BestStagingTot = BestStagingTot,
        BestStagingDown = BestStagingDown,
        BestStagingUp = BestStagingUp,
        NodesStg = InfStages,
        NodesStg_Ext = InfStages_Ext,
        NormUp = NormUp,
        NormTot = NormTot,
        NormDown = NormDown,
        CellStg = lapply(InfStages, function(x){x[as.character(ProcStruct$Partition)]}),
        CellStg_Ext = lapply(InfStages_Ext, function(x){x[as.character(ProcStruct$Partition)]})
      )
    )

  }












#
#
#
#
#   GeneExprMat.DF <- data.frame(DataStruct$ExpMat[names(ProcStruct$CellsPT), ])
#
#   DF <- PartitionData(X = DataStruct$Analysis$FinalExpMat,
#                       NodePositions = DataStruct$Analysis$FinalStruct$NodePositions,
#                       SquaredX = rowSums(DataStruct$Analysis$FinalExpMat^2))
#
#   GeneExprMat.DF.Split <- split(GeneExprMat.DF, DF$Partition)
#   GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
#   GeneExprMat.DF.Split.Sd <- lapply(GeneExprMat.DF.Split, function(x){apply(x, 2, sd)})
#
#   GeneExprMat.DF.Split.Mean.Bind <- sapply(GeneExprMat.DF.Split.Mean, cbind)
#   GeneExprMat.DF.Split.Sd.Bind <- sapply(GeneExprMat.DF.Split.Sd, cbind)
#
#   rownames(GeneExprMat.DF.Split.Mean.Bind) <- colnames(GeneExprMat.DF)
#   colnames(GeneExprMat.DF.Split.Mean.Bind) <- names(GeneExprMat.DF.Split.Mean)
#
#   rownames(GeneExprMat.DF.Split.Sd.Bind) <- colnames(GeneExprMat.DF)
#   colnames(GeneExprMat.DF.Split.Sd.Bind) <- names(GeneExprMat.DF.Split.Sd)
#
#   Reord <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
#   Reord <- Reord[Reord %in% colnames(GeneExprMat.DF.Split.Mean.Bind)]
#
#   GeneExprMat.DF.Split.Mean.Bind <- GeneExprMat.DF.Split.Mean.Bind[,Reord]
#   GeneExprMat.DF.Split.Sd.Bind <- GeneExprMat.DF.Split.Sd.Bind[,Reord]
#
#   MedianCV <- apply(GeneExprMat.DF.Split.Sd.Bind/GeneExprMat.DF.Split.Mean.Bind, 1, median, na.rm=TRUE)
#
#   # NormExp <- ProcStruct$NodesExp
#
#   NormExp <- GeneExprMat.DF.Split.Mean.Bind
#   dim(NormExp)
#
#
#
#
#
#
#
#   MedianCV <- MedianCV[rowSums(NormExp>0) > MinNodes]
#   NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
#   dim(NormExp)
#
#   SDVect <- apply(NormExp, 1, sd)
#   SDVect[is.na(SDVect)] <- 0
#   MedianCV[is.na(MedianCV)] <- 0
#
#   NormExp <- NormExp[SDVect > quantile(SDVect, QuantSel) & MedianCV < MedianCVThr,]
#
#   NormExp[NormExp <= FiltMax] <- 0
#   NormExp <- NormExp[rowSums(NormExp>0) > MinNodes, ]
#
#   NormExp <- (NormExp - apply(NormExp, 1, min))/(apply(NormExp, 1, max) - apply(NormExp, 1, min))
#
#   pheatmap::pheatmap(NormExp, cluster_cols = FALSE)
#
#   FiltStagesOnNodes <- list()
#
#
#   AVGenAct_NoNorm <- apply(NormExp[rownames(NormExp) %in% CCGenes,], 2, mean)
#   AVGenAct_NoNorm.Comp <- apply(NormExp[!(rownames(NormExp) %in% CCGenes),], 2, mean)
#
#   NormExp <- (NormExp - apply(NormExp, 1, min))/(apply(NormExp, 1, max) - apply(NormExp, 1, min))
#
#
#
#
#
#
#
#
#
#   temp <- lapply(1:length(StageInfo.UP), function(i){
#
#     if(is.null(StageInfo.UP[[i]])){
#       return(NULL)
#     }
#
#     if(sum(rownames(NormExp) %in% StageInfo.UP[[i]]) == 0){
#       return(NULL)
#     }
#
#     pheatmap::pheatmap(1*NormExp[rownames(NormExp) %in% StageInfo.UP[[i]],],
#                        show_colnames = TRUE,
#                        show_rownames = FALSE,
#                        cluster_cols = FALSE,
#                        cluster_rows = TRUE,
#                        main = paste(names(StageInfo.UP)[i], "UP"))
#   })
#
#
#
#
#
#
#
#
#   temp <- lapply(1:length(StageInfo.DOWN), function(i){
#
#     if(is.null(StageInfo.DOWN[[i]])){
#       return(NULL)
#     }
#
#     if(sum(rownames(NormExp) %in% StageInfo.DOWN[[i]]) == 0){
#       return(NULL)
#     }
#
#     pheatmap::pheatmap(1*NormExp[rownames(NormExp) %in% StageInfo.DOWN[[i]],],
#                        show_colnames = TRUE,
#                        show_rownames = FALSE,
#                        cluster_cols = FALSE,
#                        cluster_rows = TRUE,
#                        main = paste(names(StageInfo.DOWN)[i], "DOWN"))
#   })
#
#
#
#
#
#   # PercOnNodes.UP <- sapply(StageInfo.UP, function(x){
#   #   # apply(NormExp.Mod[rownames(NormExp.Mod) %in% intersect(x, Selected),], 2, sum)
#   #   # apply(NormExp[rownames(NormExp) %in% x,], 2, sum)
#   #   tMat <- NormExp[rownames(NormExp) %in% x,]
#   #   tMat[tMat < Thr] <- 0
#   #   apply(tMat, 2, mean)
#   # })
#   #
#   #
#   # PercOnNodes.DOWN <- sapply(StageInfo.DOWN, function(x){
#   #   # apply(NormExp.Mod[rownames(NormExp.Mod) %in% intersect(x, Selected),], 2, sum)
#   #   # apply(NormExp[rownames(NormExp) %in% x,], 2, sum)
#   #   tMat <- 1 - NormExp[rownames(NormExp) %in% x,]
#   #   tMat[tMat < Thr] <- 0
#   #   apply(tMat, 2, mean)
#   # })
#
#
#   PercOnNodes <- sapply(1:length(StageInfo.UP), function(i){
#
#     if(sum(rownames(NormExp) %in% union(StageInfo.UP[[i]], StageInfo.DOWN[[i]])) <= 1){
#       return(rep(0, ncol(NormExp)))
#     }
#
#     tMat.UP <- NormExp[rownames(NormExp) %in% StageInfo.UP[[i]],]
#     tMat.UP[tMat.UP < Thr] <- 0
#
#     tMat.DOWN <- 1 - NormExp[rownames(NormExp) %in% StageInfo.DOWN[[i]],]
#     tMat.DOWN[tMat.DOWN < Thr] <- 0
#
#     CombMat <- rbind(tMat.DOWN, tMat.UP)
#     colnames(CombMat) <- colnames(tMat.UP)
#
#     apply(CombMat, 2, mean)
#   })
#
#   colnames(PercOnNodes) <- names(StageInfo.UP)
#
#   barplot(t(PercOnNodes), beside = TRUE)
#
#   PercOnNodes <- PercOnNodes[, !apply(PercOnNodes == 0, 2, all)]
#
#   if(ComputeG0){
#
#     AVGenAct <- apply(NormExp[rownames(NormExp) %in% CCGenes,], 2, mean)
#     # AVGenAct.Comp <- apply(NormExp[!(rownames(NormExp) %in% CCGenes),], 2, mean)
#
#     barplot(rbind(AVGenAct_NoNorm, AVGenAct_NoNorm.Comp), beside = TRUE)
#     abline(h = G0Level)
#
#     ToKeep <- AVGenAct_NoNorm > G0Level
#
#     PercOnNodes.Comb <- t(PercOnNodes)
#     PercOnNodes.Comb[] <- 0
#     PercOnNodes.Comb[,ToKeep] <- t(PercOnNodes[ToKeep, ])/apply(PercOnNodes[ToKeep, ], 2, quantile, .9)
#
#     # G0Val <- min(apply(PercOnNodes.Comb[,ToKeep], 2, max))
#     # PercOnNodes.Comb <- rbind(G0Val*(!ToKeep), PercOnNodes.Comb)
#
#     # PercOnNodes.Comb <- rbind((1 - AVGenAct)*(!ToKeep), PercOnNodes.Comb)
#
#     PercOnNodes.Comb <- rbind(
#       (1 - AVGenAct)/quantile(1 - AVGenAct, .9),
#       PercOnNodes.Comb
#     )
#
#     PercOnNodes.Comb[1, ToKeep] <- 0
#     rownames(PercOnNodes.Comb)[1] <- "G0"
#
#     barplot(PercOnNodes.Comb, beside = TRUE)
#
#     OutPos <- apply(PercOnNodes.Comb[, ToKeep], 1, scater::isOutlier) %>%
#       apply(., 1, any)
#
#     OutVect <- ToKeep
#     OutVect[OutVect] <- OutPos
#
#     NonOutMax <- apply(PercOnNodes.Comb[-1,!OutVect & ToKeep], 1, max)
#     NonOutMax[NonOutMax == 0] <- 1
#
#     NormVect <- sapply(1:length(OutVect), function(i){
#       if(!OutVect[i]){
#         return(1)
#       } else {
#         return(max(PercOnNodes.Comb[-1, i]/NonOutMax))
#       }
#     })
#
#     PercOnNodes.Comb <- t(t(PercOnNodes.Comb)/NormVect)
#
#     barplot(PercOnNodes.Comb, beside = TRUE)
#
#   } else {
#
#     PercOnNodes.Comb <- t(PercOnNodes)/apply(PercOnNodes, 2, quantile, .9)
#
#     barplot(PercOnNodes.Comb, beside = TRUE)
#
#   }
#
#
#
#   # G0Perc <- colSums(!PlotMat[rownames(PlotMat) %in% unlist(StageInfo, use.names = FALSE),])/sum(rownames(PlotMat) %in% unlist(StageInfo, use.names = FALSE))
#
#   # PercOnNodes.UP <- t(PercOnNodes.UP)/apply(PercOnNodes.UP, 2, quantile, .9)
#   # PercOnNodes.DOWN <- t(PercOnNodes.DOWN)/apply(PercOnNodes.DOWN, 2, quantile, .9)
#   #
#   # PercOnNodes.UP[is.na(PercOnNodes.UP)] <- 0
#   # PercOnNodes.DOWN[is.na(PercOnNodes.DOWN)] <- 0
#
#   # barplot(PercOnNodes, beside = TRUE)
#
#
#
#
#   # PercOnNodes.Comb <- t(PercOnNodes)
#
#
#   # PercOnNodes <- PercOnNodes + min(PercOnNodes)
#
#   # PercOnNodes[PercOnNodes > 1] <- 1
#
#
#
#   # LowHigh <- apply(PercOnNodes, 1, quantile, .7)
#   #
#   # # barplot(PercOnNodes, beside = TRUE)
#   #
#   # pheatmap::pheatmap((PercOnNodes > LowHigh)*1, cluster_cols = FALSE, cluster_rows = FALSE)
#   #
#   # lapply(1:nrow(PercOnNodes), function(i){
#   #
#   #   PercOnNodes[i,] > LowHigh[1,i]
#   #
#   # })
#
#   pheatmap::pheatmap(PercOnNodes.Comb,
#                      cluster_cols = FALSE,
#                      cluster_rows = FALSE)
#
#   apply(PercOnNodes.Comb, 2, which.max)
#
#   # PercOnNodes[PercOnNodes < .02] <- 0
#   #
#   # pheatmap::pheatmap(PercOnNodes,
#   #                    cluster_cols = FALSE,
#   #                    cluster_rows = FALSE)
#   #
#   # apply(PercOnNodes, 2, which.max)
#
#   BestStaging <- AssociteNodes(PerMat = PercOnNodes.Comb, Mode = Mode)
#
#   TB <- apply(BestStaging, 2, function(x){table(factor(x, levels = 1:nrow(PercOnNodes.Comb)))})
#
#   rownames(TB) <- rownames(PercOnNodes.Comb)
#
#   pheatmap::pheatmap(TB, cluster_cols = FALSE, cluster_rows = FALSE,
#                      main = Title)
#
#   BestTB <- TB
#
#   InferredStages <- rownames(PercOnNodes.Comb)[apply(TB, 2, which.max)]
#   names(InferredStages) <- colnames(PercOnNodes.Comb)
#
#   CellStages <- InferredStages[paste(TaxVect)]
#   names(CellStages) <- names(TaxVect)
#
#   if(ComputeG0){
#     CellStages <- factor(CellStages, levels = c('G0', names(StageInfo.UP)))
#   } else {
#     CellStages <- factor(CellStages, levels = names(StageInfo.UP))
#   }
#
#
#   if(!is.null(names(DataStruct$Cats))){
#     TB1 <- table(DataStruct$Cats[names(CellStages)], CellStages)
#   } else {
#     TB1 <- table(DataStruct$Cats, CellStages)
#   }
#
#   print(TB1/rowSums(TB1))
#
#   pheatmap::pheatmap(TB1/rowSums(TB1), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
#                      color = rev(heat.colors(25)))
#
#   if(nrow(TB1) > 1){
#
#     print(t(t(TB1)/colSums(TB1)))
#     pheatmap::pheatmap(t(t(TB1)/colSums(TB1)), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
#                        color = rev(heat.colors(25)))
#
#   }
#
#
#   ExtStages <- c(InferredStages[length(InferredStages)], InferredStages, InferredStages[1])
#   CombStages <- rep(NA, length(InferredStages))
#   names(CombStages) <- names(InferredStages)
#
#   for(i in 2:(length(ExtStages)-1)){
#
#     if(ExtStages[i-1] == ExtStages[i] &
#        ExtStages[i+1] == ExtStages[i]){
#       CombStages[i-1] <- ExtStages[i]
#       next()
#     }
#
#     if(ExtStages[i-1] == ExtStages[i] &
#        ExtStages[i+1] != ExtStages[i]){
#       CombStages[i-1] <- paste(ExtStages[i], ExtStages[i+1], sep = " / ")
#       next()
#     }
#
#     if(ExtStages[i-1] != ExtStages[i] &
#        ExtStages[i+1] == ExtStages[i]){
#       CombStages[i-1] <- paste(ExtStages[i-1], ExtStages[i], sep = " / ")
#       next()
#     }
#
#     if(ExtStages[i-1] != ExtStages[i] &
#        ExtStages[i+1] != ExtStages[i]){
#       CombStages[i-1] <- paste(ExtStages[i-1], ExtStages[i], ExtStages[i+1], sep = " / ")
#       next()
#     }
#
#   }
#
#
#   CellStages_Ext <- CombStages[as.character(TaxVect)]
#   names(CellStages_Ext) <- names(TaxVect)
#   # CellStages_Ext <- factor(CellStages_Ext)
#
#
#   if(!is.null(names(DataStruct$Cats))){
#     TB2 <- table(DataStruct$Cats[names(CellStages_Ext)], CellStages_Ext)
#   } else {
#     TB2 <- table(DataStruct$Cats, CellStages_Ext)
#   }
#
#   print(TB2/rowSums(TB2))
#
#   pheatmap::pheatmap(TB2/rowSums(TB2), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
#                      color = rev(heat.colors(25)))
#
#   if(nrow(TB2) > 1){
#
#     print(t(t(TB2)/colSums(TB2)))
#     pheatmap::pheatmap(t(t(TB2)/colSums(TB2)), cluster_rows = FALSE, cluster_cols = FALSE, main = Title,
#                        color = rev(heat.colors(25)))
#
#
#   }
#
#
#   return(list(Inferred = InferredStages,
#               Extinferred = CombStages,
#               CellStages = CellStages,
#               CellStages_Ext = CellStages_Ext,
#               StageMat = PercOnNodes.Comb,
#               TB1 = TB1, TB2 = TB2,
#               BestTB = BestTB))

}































#' Title
#'
#' @param StageMatrix
#' @param NodePenalty
#'
#' @return
#' @export
#'
#' @examples
FitStagesCirc <- function(StageMatrix, NodePenalty, Mode = 1) {

  NormStageMatrix <- StageMatrix

  # Find which columns are associated with at least one stage

  ToAnalyze <- which(colSums(StageMatrix)>0)

  if(length(ToAnalyze)<nrow(StageMatrix)){

    # Not enough columns. Padding around them

    PaddedIdxs <- c(
      (min(ToAnalyze) - nrow(StageMatrix)):(min(ToAnalyze)-1),
      ToAnalyze,
      (max(ToAnalyze)+1):(max(ToAnalyze) + nrow(StageMatrix))
    )

    PaddedIdxs[PaddedIdxs > ncol(StageMatrix)] <-  PaddedIdxs[PaddedIdxs > ncol(StageMatrix)] - ncol(StageMatrix)
    PaddedIdxs[PaddedIdxs <= 0] <-  PaddedIdxs[PaddedIdxs <= 0] + ncol(StageMatrix)

    ToAnalyze <- sort(unique(PaddedIdxs))

  }

  # Selecting the matrix that will be analysed and saving the indices

  NormStageMatrix <- NormStageMatrix[,ToAnalyze]
  NormStageMatrixIdx <- ToAnalyze

  NormNodePenalty <- NodePenalty[ToAnalyze]

  Possibilities <- combn(c(1:ncol(NormStageMatrix), rep(NA, nrow(NormStageMatrix))), nrow(NormStageMatrix))

  NoChange <- apply(is.na(Possibilities), 2, sum) == nrow(NormStageMatrix)

  ToKeep <- (!apply(Possibilities, 2, is.unsorted, na.rm = TRUE))
  #   (apply(!is.na(Possibilities), 2, sum) != 1) &
  #   (!is.na(Possibilities[nrow(NormStageMatrix),]))

  Possibilities <- Possibilities[, ToKeep | NoChange]
  dim(Possibilities) <- c(length(Possibilities)/(sum(ToKeep | NoChange)),
                          sum(ToKeep | NoChange))

  PathPenality <- function(ChangeNodes, InitialStage, Mode) {

    Sphases <- rep(InitialStage, ncol(NormStageMatrix))

    tStart <- NA
    tEnd <- NA

    for (i in 1:(length(ChangeNodes)-1)) {

      if(!is.na(ChangeNodes[i])){
        tStart <- ChangeNodes[i]
      }

      tEnd <- ChangeNodes[i+1]

      if(is.na(tStart) | is.na(tEnd)){
        next()
      }

      Sphases[(tStart+1):(tEnd)] <- InitialStage + i
    }

    Sphases[Sphases>length(ChangeNodes)] <- Sphases[Sphases>length(ChangeNodes)] - length(ChangeNodes)

    if(Mode == 1){
      # Squared distance from the maximum
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, sum) - diag(NormStageMatrix[Sphases,]))
        )
      )
    }

    if(Mode == 2){
      # Sum of "off-stage" contributions
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, max) - mapply("[[", apply(NormStageMatrix, 2, as.list), Sphases))
        )
      )
    }

    if(Mode == 3){
      # Sum binary difference from maximum
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, which.max) != Sphases)
        )
      )
    }

    if(Mode == 4){
      # Sum binary difference from maximum
      tDiff <- abs(apply(NormStageMatrix, 2, which.max) - Sphases)
      tDiff[tDiff > nrow(NormStageMatrix)/2] <- nrow(NormStageMatrix) - tDiff[tDiff > nrow(NormStageMatrix)/2]
      return(
        sum(
          NormNodePenalty*tDiff
        )
      )
    }

  }

  CombinedInfo <- NULL

  for (i in 1:nrow(NormStageMatrix)) {
    CombinedInfo <- cbind(CombinedInfo,
                          rbind(rep(i, sum(ToKeep | NoChange)),
                                apply(Possibilities, 2, PathPenality, InitialStage = i, Mode = Mode),
                                1:(sum(ToKeep | NoChange))
                          )
    )
  }

  return(list(Penality = CombinedInfo, Possibilities = Possibilities))

}





















#' Title
#'
#' @param PerMat
#'
#' @return
#' @export
#'
#' @examples
AssociteNodes <- function(PerMat, Mode = 1) {

  tictoc::tic()
  print("Direct staging")
  Staging <- FitStagesCirc(StageMatrix = PerMat,
                           NodePenalty = rep(1, ncol(PerMat)),
                           Mode = Mode)
  tictoc::toc()

  tictoc::tic()
  print("Reverse staging")
  StagingRev <- FitStagesCirc(StageMatrix = PerMat[, rev(1:ncol(PerMat))],
                              NodePenalty = rep(1, ncol(PerMat)),
                              Mode = Mode)
  tictoc::toc()

  AllPenality <- rbind(cbind(Staging$Penality, StagingRev$Penality), rep(1:2, each=ncol(Staging$Penality)))

  Idxs <- which(AllPenality[2, ] == min(AllPenality[2, ]))

  SelPenality <- AllPenality[,Idxs]
  dim(SelPenality) <- c(4, length(SelPenality)/4)

  DirectPenality <- NULL
  DirectChanges <- NULL

  if(sum(SelPenality[4,] == 1) > 0){

    ExpandStages <- function(idx) {

      ChangeNodes <- Staging$Possibilities[ , SelPenality[3, idx]]

      StageVect <- rep(SelPenality[1, idx], ncol(PerMat))

      tStart <- NA
      tEnd <- NA

      for (i in 1:(length(ChangeNodes)-1)) {

        if(!is.na(ChangeNodes[i])){
          tStart <- ChangeNodes[i]
        }

        tEnd <- ChangeNodes[i+1]

        if(is.na(tStart) | is.na(tEnd)){
          next()
        }

        StageVect[(tStart+1):(tEnd)] <- SelPenality[1, idx] + i
      }

      StageVect[StageVect>length(ChangeNodes)] <- StageVect[StageVect>length(ChangeNodes)] - length(ChangeNodes)

      return(StageVect)

    }

    SelPenIdx <- which(SelPenality[4,]==1)

    DirectPenality <- SelPenality[2,SelPenIdx]
    DirectChanges <- t(sapply(SelPenIdx, ExpandStages))

  }

  ReversePenality <- NULL
  ReverseChanges <- NULL

  if(sum(SelPenality[4,] == 2) > 0){

    ExpandStages <- function(idx) {

      ChangeNodes <- StagingRev$Possibilities[ , SelPenality[3, idx]]

      StageVect <- rep(SelPenality[1, idx], ncol(PerMat))

      tStart <- NA
      tEnd <- NA

      for (i in 1:(length(ChangeNodes)-1)) {

        if(!is.na(ChangeNodes[i])){
          tStart <- ChangeNodes[i]
        }

        tEnd <- ChangeNodes[i+1]

        if(is.na(tStart) | is.na(tEnd)){
          next()
        }

        StageVect[(tStart+1):(tEnd)] <- SelPenality[1, idx] + i
      }

      StageVect[StageVect>length(ChangeNodes)] <- StageVect[StageVect>length(ChangeNodes)] - length(ChangeNodes)

      return(rev(StageVect))

    }

    SelPenIdx <- which(SelPenality[4,]==2)

    ReversePenality <- SelPenality[2,SelPenIdx]
    ReverseChanges <- t(sapply(SelPenIdx, ExpandStages))

  }

  AllStg <- rbind(DirectChanges, ReverseChanges)
  AllPen <- c(DirectPenality, ReversePenality)
  AllDir <- c(rep("Dir", length(DirectPenality)),
              rep("Rev", length(ReversePenality)))
  colnames(AllStg) <- colnames(PerMat)

  return(AllStg)

}







