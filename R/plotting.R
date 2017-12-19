#################################################################################
# PlotOnStages ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param Structure
#' @param TaxonList
#' @param Categories
#' @param PrinGraph
#' @param Net
#' @param SelThr
#' @param ComputeOverlaps
#' @param RotatioMatrix
#' @param ExpData
#'
#' @return
#' @export
#'
#' @examples
PlotOnStages <- function(Structure,
                         PrinStruct,
                         ExpData,
                         SelThr = NULL,
                         SmoothPoints.Growth = 1,
                         SmoothPoints.Shrink = 1,
                         ComputeOverlaps = TRUE,
                         nGenes = 10,
                         OrderOnCat = FALSE,
                         MinCellPerNode = 3,
                         Title = '') {

  # Structure <- "Circle"
  # PrinStruct <- Data.Kowa.CV$Analysis
  # ExpData <- Data.Kowa.CV$ExpMat
  # SelThr = .3
  # SmoothPoints.Growth = 1
  # SmoothPoints.Shrink = 1
  # OrderOnCat = TRUE

  # Check categories
  Categories <- PrinStruct$FinalGroup

  if(!is.factor(Categories)){
    stop("Categories must be a factor")
  }

  if(is.null(SelThr)){
    SelThr <- 1.1*(1/length(levels(Categories)))
  }

  # Get Info on the nodes

  Nodes <- PrinStruct$FinalStruct$NodePositions
  nNodes <- nrow(Nodes)

  # Construct graph
  Net <- ElPiGraph.R::ConstructGraph(PrintGraph = PrinStruct$FinalStruct)

  # Find the best Staging if circular

  Partition <- ElPiGraph.R::PartitionData(X = PrinStruct$FinalExpMat, NodePositions = Nodes, SquaredX = rowSums(PrinStruct$FinalExpMat^2))$Partition
  ProjStruct <- ElPiGraph.R::project_point_onto_graph(X = PrinStruct$FinalExpMat, NodePositions = Nodes,
                                                   Edges = PrinStruct$FinalStruct$Edges$Edges, Partition = Partition)

  TB <- table(Categories, factor(Partition, levels = paste(1:nNodes)))
  # colnames(TB)

  print("Step I - Finding the best path")

  StepIDone <- FALSE

  if(Structure == 'Circle'){

    StepIDone <- TRUE

    print("Getting all circular subisomorphisms")

    AllPaths <- igraph::subgraph_isomorphisms(pattern = igraph::graph.lattice(nNodes, directed = FALSE, circular = TRUE), target = Net)

    if(OrderOnCat){

      SummInfo <- sapply(AllPaths, function(subIso){

        Pt <- ElPiGraph.R::getPseudotime(Edges = PrinStruct$FinalStruct$Edges$Edges, ProjStruct = ProjStruct,
                                EdgeSeq = c(names(subIso), names(subIso)[1]))

        AGG <- aggregate(Pt$Pt, by = list(Categories), median)
        AGG2 <- aggregate(Pt$Pt, by = list(Categories), min)
        AGG3 <- aggregate(Pt$Pt, by = list(Categories), max)

        Sorted <- FALSE

        if(S4Vectors::isSorted(AGG[,2])){
          return(c(1, summary(aov(Pt$Pt ~ Categories))[[1]][1,"Pr(>F)"], AGG2[1,2], AGG[1,2]))
        } else {
          return(c(2, summary(aov(Pt$Pt ~ Categories))[[1]][1,"Pr(>F)"], AGG2[1,2], AGG[1,2]))
        }

      })

      print("Finding the best path")

      SummInfo <- rbind(1:ncol(SummInfo), SummInfo)

      if(any(SummInfo[2,] == 1)){
        SummInfo <- SummInfo[ , SummInfo[2,] == 1]
        print("Strongly consecutive stages found")
      } else {
        print("Unable to find strongly consecutive stages")
      }

      if(length(SummInfo) == 5){
        SummInfo <- matrix(SummInfo, ncol = 1)
      }

      Selected <- SummInfo[, SummInfo[3, ] == min(SummInfo[3, ])]

      print(Selected)

      if(length(Selected)>5){
        Selected <- Selected[,which.min(Selected[4,])]
        print(Selected)
      }

      SelPath <- AllPaths[[Selected[1]]]

      Pt <- ElPiGraph.R::getPseudotime(Edges = PrinStruct$FinalStruct$Edges$Edges, ProjStruct = ProjStruct,
                                    EdgeSeq = c(names(SelPath), names(SelPath)[1]))

      boxplot(Pt$Pt ~ Categories, main = Title, ylab = "Pseudotime on path")

    } else {
      SelPath <- AllPaths[[sample(1:length(AllPaths), 1)]]
    }

  }

  # if(Structure == 'Line'){
  #
  #   StepIDone <- TRUE
  #
  #   print("Getting all line subisomorphisms")
  #
  #   AllPaths <- GetLongestPath(Net = Net, Structure = Structure, Circular = FALSE)
  #
  #   if(OrderOnCat){
  #
  #     SummInfo <- NULL
  #
  #     for(i in 1:nrow(AllPaths$VertNumb)){
  #       SelPath <- AllPaths$VertNumb[i,]
  #       # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
  #       # SelPath <- SelPath[-length(SelPath)]
  #
  #       Reordered <- TaxVect
  #       for(j in 1:length(SelPath)){
  #         Reordered[TaxVect == SelPath[j]] <- j
  #       }
  #
  #       NumReord <- as.numeric(as.character(Reordered))
  #
  #       AGG <- aggregate(NumReord, by = list(Categories), median)
  #       AGG2 <- aggregate(NumReord, by = list(Categories), min)
  #       AGG3 <- aggregate(NumReord, by = list(Categories), max)
  #
  #       if(S4Vectors::isSorted(AGG[,2])){
  #         SummInfo <- rbind(SummInfo,
  #                           c(i, 1, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
  #                             AGG2[1,2], AGG[1,2])
  #         )
  #
  #         # boxplot(NumReord ~ Categories, main = i)
  #
  #       }
  #
  #
  #       if(S4Vectors::isSorted(rev(AGG[,2]))){
  #         SummInfo <- rbind(SummInfo,
  #                           c(i, 2, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
  #                             nNodes - AGG3[1,2] + 1,
  #                             nNodes - AGG[1,2] + 1)
  #         )
  #         # boxplot(Reordered ~ Categories, main = paste(i, "rev"))
  #       }
  #     }
  #
  #     print("Finding the best path")
  #
  #     Selected <- SummInfo[SummInfo[, 3] == min(SummInfo[, 3]), ]
  #
  #     if(length(Selected)>5){
  #       Selected <- Selected[which.min(Selected[,4]),]
  #     }
  #
  #     SelPath <- AllPaths$VertNumb[Selected[1],]
  #     # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
  #     SelPath <- SelPath[-length(SelPath)]
  #
  #     if(Selected[2] == 2){
  #       SelPath <- rev(SelPath)
  #     }
  #
  #     Reordered <- TaxVect
  #     for(j in 1:length(SelPath)){
  #       Reordered[as.character(TaxVect) == SelPath[j]] <- j
  #     }
  #
  #     NumReord <- as.numeric(as.character(Reordered))
  #
  #     boxplot(NumReord ~ Categories)
  #
  #   } else {
  #     SelPath <- AllPaths$VertNumb[sample(x = 1:nrow(AllPaths$VertNumb), size = 1),]
  #   }
  #
  #
  # }

  if(!StepIDone){
    stop("Unsupported Structure")
  }

  print("Step II")

  ExtendedTB <- TB[,names(SelPath)]
  
  ExtendedTB <- matrix(ExtendedTB, ncol = length(SelPath))
  
  colnames(ExtendedTB) <- names(SelPath)
  rownames(ExtendedTB) <- rownames(TB)

  if(OrderOnCat & nrow(ExtendedTB)>1){
    barplot(t(t(ExtendedTB)/colSums(ExtendedTB)), col = rainbow(nrow(TB)), beside = TRUE, las = 2,
            legend.text = rownames(ExtendedTB), args.legend = list(x = "top", fill=rainbow(nrow(TB))),
            ylim = c(0, 1.25), yaxt = "n")
    axis(2, seq(from=0, to=1, by=.25), las=2)
  }

  print(ExtendedTB)

  # SelPath <- SelPath[-length(SelPath)]
  SelPathSTG <- rep(NA, length(SelPath))
  
  
  PercMat <- t(t(ExtendedTB)/colSums(ExtendedTB))
  PercMat[is.na(PercMat)] <- 0
  PercMat[,colSums(ExtendedTB) < MinCellPerNode] <- 0
  BinPercMat <- (PercMat > SelThr)



  if(OrderOnCat){

    # Fill gaps with TRUE

    # if(Structure == 'Circle'){
    #
    #   SelGrop <- apply(BinPercMat, 1, function(x){
    #     range(which(x))})
    #
    #   for(j in 1:ncol(SelGrop)){
    #     if(any(is.infinite(SelGrop[,j]))){
    #       next
    #     } else {
    #       Diff <- SelGrop[1,j] + ncol(BinPercMat) + 1 - SelGrop[2,j] - 2
    #
    #       if(Diff > 0 & Diff <= SmoothPoints){
    #         BinPercMat[colnames(SelGrop)[j],1:(SelGrop[1,j]-1)] <- TRUE
    #         BinPercMat[colnames(SelGrop)[j],(SelGrop[2,j]+1):ncol(BinPercMat)] <- TRUE
    #       }
    #     }
    #   }
    #
    # }

    if(Structure == 'Circle'){
      BinPercMat <- cbind(BinPercMat, BinPercMat[,1])
    }

    WhichStage <- apply(BinPercMat, 1, which)

    WhichStage.Diff <- lapply(WhichStage, function(x){x[-1] - x[-length(x)]})

    lapply(as.list(1:length(WhichStage.Diff)), function(j){
      GapsToFill <- which(WhichStage.Diff[[j]] > 1 & WhichStage.Diff[[j]] <= SmoothPoints.Growth + 1)

      for(k in GapsToFill){
        BinPercMat[names(WhichStage)[j], WhichStage[[j]][k]:WhichStage[[j]][k+1]] <<- TRUE
      }

    })


    if(Structure == 'Circle'){
      BinPercMat <- BinPercMat[,-ncol(BinPercMat)]
    }


    # Fill gaps with FALSE

    # if(Structure == 'Circle'){
    #
    #   SelGrop <- apply(BinPercMat, 1, function(x){
    #     range(which(!x))})
    #
    #   for(j in 1:ncol(SelGrop)){
    #     if(any(is.infinite(SelGrop[,j]))){
    #       next
    #     } else {
    #       Diff <- SelGrop[1,j] + ncol(BinPercMat) + 1 - SelGrop[2,j] - 2
    #
    #       if(Diff > 0 & Diff <= SmoothPoints){
    #         BinPercMat[colnames(SelGrop)[j],1:(SelGrop[1,j]-1)] <- FALSE
    #         BinPercMat[colnames(SelGrop)[j],(SelGrop[2,j]+1):ncol(BinPercMat)] <- FALSE
    #       }
    #     }
    #   }
    #
    # }

    if(Structure == 'Circle'){
      BinPercMat <- cbind(BinPercMat, BinPercMat[,1])
    }

    WhichStage <- apply(BinPercMat, 1, function(x){which(!x)})

    WhichStage.Diff <- lapply(WhichStage, function(x){x[-1] - x[-length(x)]})

    lapply(as.list(1:length(WhichStage.Diff)), function(j){
      GapsToFill <- which(WhichStage.Diff[[j]] > 1 & WhichStage.Diff[[j]] <= SmoothPoints.Shrink + 1)

      for(k in GapsToFill){
        BinPercMat[names(WhichStage)[j], WhichStage[[j]][k]:WhichStage[[j]][k+1]] <<- FALSE
      }

    })

    if(Structure == 'Circle'){
      BinPercMat <- BinPercMat[,-ncol(BinPercMat)]
    }

    print(1*BinPercMat)

    if(ComputeOverlaps){

      # constructing stage overlaps

      OvefLapCat <- list()

      for(i in 1:nrow(BinPercMat)){

        if(i == nrow(BinPercMat)){
          if(Structure == 'Circle'){
            OvefLapCat[[i]] <- BinPercMat[i,] & BinPercMat[1,]
          }
        } else {
          OvefLapCat[[i]] <- BinPercMat[i,] & BinPercMat[i+1,]
        }

      }

      BinPercMatExt <- NULL

      for(i in 1:nrow(BinPercMat)){

        if(i == 1){
          if(Structure == 'Circle'){
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[1]] & !OvefLapCat[[nrow(BinPercMat)]],
                                   OvefLapCat[[i]])
          } else {
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[1]],
                                   OvefLapCat[[i]])
          }
          next()
        }

        if(i == nrow(BinPercMat)){
          if(Structure == 'Circle'){
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i-1]] & !OvefLapCat[[i]],
                                   OvefLapCat[[i]])
          } else {
            BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i-1]])
          }
          next()
        }

        BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i]] & !OvefLapCat[[i-1]],
                               OvefLapCat[[i]])

      }

      OverlapCat <- paste(levels(Categories)[-length(levels(Categories))],
                          levels(Categories)[-1], sep = "+")

      if(Structure == 'Circle'){
        OverlapCat <- c(OverlapCat,
                        paste(levels(Categories)[length(levels(Categories))],
                              levels(Categories)[1], sep = "+")
        )
      }

      if(Structure == 'Circle'){
        rownames(BinPercMatExt) <- as.vector(rbind(rownames(BinPercMat), OverlapCat))
      } else {
        rownames(BinPercMatExt) <- as.vector(rbind(rownames(BinPercMat), c(OverlapCat, NA)))[-nrow(BinPercMat)*2]
      }

    } else {
      BinPercMatExt <- BinPercMat
    }

    print(1*BinPercMatExt)

    AllCat <- rownames(BinPercMatExt)

    # if(Structure == "Circle"){
    #   ExtPath <- ExtPath[-length(ExtPath)]
    #   BinPercMatExt <- BinPercMatExt[, -ncol(BinPercMatExt)]
    # }

    LowStages <- 1
    if(ComputeOverlaps){
      LowStages <- 1:2
    }


    BinPercMatExt.ColNames <-  colnames(BinPercMatExt)
    while(any(BinPercMatExt[LowStages,ncol(BinPercMatExt)]) &
          all(!BinPercMatExt[setdiff(1:nrow(BinPercMatExt), LowStages),ncol(BinPercMatExt)])){
      BinPercMatExt <- cbind(BinPercMatExt[, ncol(BinPercMatExt)],
                             BinPercMatExt[, -ncol(BinPercMatExt)])
      BinPercMatExt.ColNames <- c(BinPercMatExt.ColNames[length(BinPercMatExt.ColNames)],
                                  BinPercMatExt.ColNames[-length(BinPercMatExt.ColNames)])
    }

    colnames(BinPercMatExt) <- BinPercMatExt.ColNames

    HighStages <- nrow(BinPercMatExt)
    if(ComputeOverlaps & Structure == 'Circle'){
      HighStages <- c(nrow(BinPercMatExt)-1, nrow(BinPercMatExt))
    }

    while(any(BinPercMatExt[HighStages,1]) &
          all(!BinPercMatExt[setdiff(1:nrow(BinPercMatExt), HighStages),1])){
      BinPercMatExt <- cbind(BinPercMatExt[, -1], BinPercMatExt[, 1])
      BinPercMatExt.ColNames <- c(BinPercMatExt.ColNames[-1],
                                  BinPercMatExt.ColNames[1])
    }

    colnames(BinPercMatExt) <- BinPercMatExt.ColNames

    # if(Structure == 'Circle'){
    #   ExtPath <- c(ExtPath, ExtPath[1])
    #   BinPercMatExt <- cbind(BinPercMatExt,
    #                          rep(FALSE, nrow(BinPercMatExt)))
    # }

  } else {
    BinPercMatExt <- BinPercMat
    AllCat <- rownames(BinPercMatExt)
  }

  if(nrow(BinPercMatExt)>1){
    Idxs <- apply(BinPercMatExt, 1, which)
  } else {
    Idxs <- list(Cat = which(BinPercMatExt))
  }

  print("Step III")
  print(Idxs)

  Bond <- lapply(Idxs[lapply(Idxs, length) > 1], range)

  if(!S4Vectors::isSorted(unlist(Bond))){
    warning("Stages are not sequential!")
  }

  print(1*BinPercMat)

  print(unlist(Bond))

  if(any(PrinStruct$PCACenter != FALSE)){
    NodeOnGenes <- t(Nodes %*% t(PrinStruct$PCARotation)) + PrinStruct$PCACenter
  } else {
    NodeOnGenes <- t(Nodes %*% t(PrinStruct$PCARotation))
  }

  # OrderedPoints <- OrderOnPath(PrinGraph = PrinGraph, Path = as.numeric(ExtPath),
  #                              PointProjections = PointProjections)


  if(nGenes > 0){
    SelIdxs <- order(apply(NodeOnGenes, 2, var), decreasing = TRUE)[1:nGenes]
  } else {
    SelIdxs <- 1:2
  }

  if(Structure == 'Circle'){
    EdgSeq <- colnames(BinPercMatExt)[c(1:ncol(BinPercMatExt), 1)]
    BinPercMatExt <- BinPercMatExt[, c(1:ncol(BinPercMatExt), 1)]
    dim(BinPercMatExt) <- c(length(BinPercMatExt)/length(EdgSeq), length(EdgSeq))
    colnames(BinPercMatExt) <- EdgSeq
    
    Pt <- ElPiGraph.R::getPseudotime(Edges = PrinStruct$FinalStruct$Edges$Edges, ProjStruct = ProjStruct,
                                  EdgeSeq = EdgSeq)
  }



  for(Idx in SelIdxs){

    p <- ggplot2::ggplot(data = data.frame(x=Pt$NodePos,
                                           y=NodeOnGenes[Idx,as.numeric(colnames(BinPercMatExt))]),
                         mapping = ggplot2::aes(x = x, y = y, color="PC")) +
      ggplot2::labs(x = "Pseudotime", y="Gene expression",
                    title = paste(Title, rownames(NodeOnGenes)[Idx]), sep=' / ') +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    RecCoord <- NULL

    if(OrderOnCat){

      for(i in 1:length(Idxs)){
        LowIds <- Idxs[[i]]-1
        LowIds[LowIds == 0] <- NA

        HiIds <- Idxs[[i]]+1
        HiIds[HiIds == ncol(BinPercMatExt) + 1] <- NA

        LowCoord <- Pt$NodePos[LowIds]
        MidCoord <- Pt$NodePos[Idxs[[i]]]
        HighCoord <- Pt$NodePos[HiIds]

        RecCoord <- rbind(RecCoord, cbind(
          colMeans(rbind(LowCoord, MidCoord)),
          MidCoord,
          colMeans(rbind(MidCoord, HighCoord)),
          rep(names(Idxs)[i], length(MidCoord))
        )
        )

      }

      colnames(RecCoord) <- c("Min", "Med", "Max", "Stage")
      RecCoord <- data.frame(RecCoord)
      RecCoord$Min <- as.numeric(as.character(RecCoord$Min))
      RecCoord$Med <- as.numeric(as.character(RecCoord$Med))
      RecCoord$Max <- as.numeric(as.character(RecCoord$Max))


      RecCoord$Stage <- factor(as.character(RecCoord$Stage), levels = AllCat)

      RecCoord$Min[is.na(RecCoord$Min)] <- RecCoord$Med[is.na(RecCoord$Min)]
      RecCoord$Max[is.na(RecCoord$Max)] <- RecCoord$Med[is.na(RecCoord$Max)]

      p <- p + ggplot2::geom_rect(data = RecCoord, mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                                  ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
        ggplot2::geom_point(data = data.frame(x=Pt$Pt,
                                              y=ExpData[rownames(PrinStruct$FinalExpMat),rownames(NodeOnGenes)[Idx]]),
                            mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))

      print(p)

    } else {

      p <- p + ggplot2::geom_point(data = data.frame(x=Pt$Pt,
                                                     y=ExpData[rownames(PrinStruct$FinalExpMat),rownames(NodeOnGenes)[Idx]]),
                            mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))

    }



  }

  CellsPT <- Pt$Pt
  names(CellsPT) <- rownames(PrinStruct$FinalExpMat)
  ExtPath <- colnames(BinPercMatExt)

  return(list(Structure = Structure,
              ExtPath = ExtPath,
              NodesPT = Pt$NodePos,
              NodesExp = NodeOnGenes[order(apply(NodeOnGenes, 1, var), decreasing = TRUE), as.numeric(ExtPath)],
              CellsPT = CellsPT,
              Partition = Partition,
              BinPercMatExt = BinPercMatExt,
              BinPercMat = BinPercMat,
              CellExp = ExpData[, order(apply(ExpData, 2, var), decreasing = TRUE)],
              RecCoord = RecCoord,
              StageOnNodes = Idxs
              )
         )

}












#' Title
#'
#' @param WorkStruct
#' @param Expression
#' @param Name
#' @param gName
#' @param SpanVal
#' @param CatOrder
#'
#' @return
#' @export
#'
#' @examples
PlotOnPseudotime <- function(WorkStruct,
                             Expression,
                             Name = '',
                             gName,
                             SpanVal=.3,
                             CatOrder = NULL,
                             legend.position = "bottom") {

  # WorkStruct = InputList[[i]]$OrderedData
  # Expression = InputList[[i]]$Expression
  # Name = InputList[[i]]$Name
  # gName = "Cdk1"
  # SpanVal = .1
  # CatOrder = NULL

  ReOrd.Sel <- match(rownames(WorkStruct$CellExp), names(WorkStruct$CellsPT))
  ReOrd.Sel.Mat <- match(rownames(Expression), names(WorkStruct$CellsPT))

  if(!is.null(CatOrder)){
    WorkStruct$RecCoord$Stage <- factor(as.character(WorkStruct$RecCoord$Stage), levels = CatOrder)
  }

  if(gName %in% rownames(WorkStruct$NodesExp)){

    SmoothData <- data.frame(x=c(WorkStruct$CellsPT[ReOrd.Sel],
                                 WorkStruct$CellsPT[ReOrd.Sel] + max(WorkStruct$NodesPT),
                                 WorkStruct$CellsPT[ReOrd.Sel] - max(WorkStruct$NodesPT)),
                             y=c(WorkStruct$CellExp[,gName],
                                 WorkStruct$CellExp[,gName],
                                 WorkStruct$CellExp[,gName]))

    p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$NodesPT,
                                            y=WorkStruct$NodesExp[gName,]),
                          mapping = ggplot2::aes(x = x, y = y, color="EpG")) +
      ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(Name, "/", gName)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = legend.position) +
      ggplot2::geom_rect(data = WorkStruct$RecCoord,
                         mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                         ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
      ggplot2::geom_smooth(data = SmoothData,
                           mapping = ggplot2::aes(x=x, y=y, color="Data"),
                           inherit.aes = FALSE, span = SpanVal, method = "loess") +
      ggplot2::geom_point(data = data.frame(x=as.vector(WorkStruct$CellsPT[ReOrd.Sel]),
                                            y=as.vector(WorkStruct$CellExp[,gName])),
                          mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
      ggplot2::geom_point() + ggplot2::geom_line() +
      ggplot2::scale_color_manual(name = "", values = c(Data = "blue", EpG = "black")) +
      ggplot2::coord_cartesian(xlim = c(0, max(WorkStruct$NodesPT)))

  } else {

    if(gName %in% colnames(Expression)){

      SmoothData <- data.frame(x=c(WorkStruct$CellsPT[ReOrd.Sel.Mat],
                                   WorkStruct$CellsPT[ReOrd.Sel.Mat] + max(WorkStruct$NodesPT),
                                   WorkStruct$CellsPT[ReOrd.Sel.Mat] - max(WorkStruct$NodesPT)),
                               y=c(unlist(Expression[,gName]),
                                   unlist(Expression[,gName]),
                                   unlist(Expression[,gName])))

      p1 <- ggplot2::ggplot(data = data.frame(x=as.vector(WorkStruct$CellsPT[ReOrd.Sel.Mat]),
                                              y=unlist(Expression[,gName])),
                            mapping = ggplot2::aes(x=x, y=y, color="Data")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = legend.position) +
        ggplot2::geom_rect(data = WorkStruct$RecCoord,
                           mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                           ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
        ggplot2::geom_smooth(data = SmoothData,
                             mapping = ggplot2::aes(x=x, y=y, color="Data"),
                             inherit.aes = FALSE, span = SpanVal, method = "loess") +
        ggplot2::geom_point(alpha=.5) +
        ggplot2::scale_color_manual(name = "", values = c(Data = "blue", EpG = "black")) +
        ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(Name, "/", gName)) +
        ggplot2::coord_cartesian(xlim = c(0, max(WorkStruct$NodesPT)))

    } else {

      p1 <- ggplot2::ggplot(data = data.frame(x=as.vector(WorkStruct$CellsPT[ReOrd.Sel.Mat]),
                                              y=unlist(Expression[,1])),
                            mapping = ggplot2::aes(x=x, y=y, color="Data")) +
        ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(Name, "/", gName)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = legend.position) +
        ggplot2::geom_rect(data = WorkStruct$RecCoord,
                           mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                           ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
        ggplot2::scale_color_manual(name = "", values = c(Data = "blue", EpG = "black")) +
        ggplot2::coord_cartesian(xlim = c(0, max(WorkStruct$NodesPT)))

    }

  }

  return(p1)

}













