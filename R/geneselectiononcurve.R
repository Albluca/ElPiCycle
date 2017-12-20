#' Title
#'
#' @param TaxVect
#' @param ExpMat
#' @param Thr
#' @param MadsThr
#' @param Mode
#' @param Net
#' @param PriGraph
#' @param Proj
#'
#' @return
#' @export
#'
#' @examples
SelectGenes <- function(Partition,
                        ExpMat,
                        Mode = "CV",
                        AggFun = min,
                        Span = .75,
                        Net = NULL,
                        Edges = NULL,
                        ProjStruct = NULL,
                        nCores = 1) {



  if(Mode == "CV" | is.null(Net) | is.null(Edges)){

    print(paste("Computing coefficient of variation on", ncol(ExpMat), "genes and", length(unique(Partition)), "pseudotime points on a single processor."))

    GeneExprMat.DF.Split <- split(data.frame(ExpMat), Partition)

    GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
    GeneExprMat.DF.Split.Sd <- lapply(GeneExprMat.DF.Split, function(x){apply(x, 2, sd)})

    GeneExprMat.DF.MeanRemoved <-
      lapply(as.list(1:length(GeneExprMat.DF.Split)), function(i){
        GeneExprMat.DF.Split.Sd[[i]]/GeneExprMat.DF.Split.Mean[[i]]
      })

    GeneExprMat.DF.MeanRemoved.All <- do.call(rbind, GeneExprMat.DF.MeanRemoved)
    colnames(GeneExprMat.DF.MeanRemoved.All) <- colnames(ExpMat)

    RetVal <- apply(abs(GeneExprMat.DF.MeanRemoved.All), 2, AggFun, na.rm = TRUE)

    return(RetVal)

  }

  if(Mode == "SmoothOnCircleNodes"){

    NodeOrder <- igraph::subgraph_isomorphisms(target = Net, pattern = igraph::make_lattice(igraph::vcount(Net), circular = TRUE))[[1]] %>%
      names(.)

    Pt <- ElPiGraph.R::getPseudotime(Edges = Edges, ProjStruct = ProjStruct, EdgeSeq = c(NodeOrder, NodeOrder[1]))

    NodesPT <- Pt$NodePos
    CellPT <- Pt$Pt
    names(CellPT) <- rownames(ExpMat)
    
    TopPT <- sort(CellPT, decreasing = FALSE)[1:ceiling(length(CellPT)*Span)]
    LowPT <- sort(CellPT, decreasing = TRUE)[1:ceiling(length(CellPT)*Span)]
    
    NewSpan <- Span/(1+2*Span)
    

    ExtPT <- c(
      TopPT + Pt$PathLen,
      CellPT,
      LowPT - Pt$PathLen
    )
    
    Sorted <- sort(ExtPT)
  
    FitFun <- function(x){
      LOE <- loess(x ~ Sorted, span = NewSpan)
      predict(LOE, data.frame(Sorted = NodesPT), se = TRUE)
    }
    
    if(nCores <= 1){
      
      print(paste("Computing loess smoothers on", ncol(ExpMat), "genes and", length(Sorted), "pseudotime points on a single processor. This may take a while ..."))
      
      tictoc::tic()
      
      op <- pboptions(type = "timer")
      PredMat <- pbapply::pbapply(ExpMat[names(Sorted),], 2, FitFun)
      pboptions(op)
      
      
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
      
      print(paste("Computing loess smoothers on", ncol(ExpMat), "genes and", length(Sorted), "pseudotime points using", nCores, "processors. This may take a while ..."))
      
      tictoc::tic()
      cl <- parallel::makeCluster(nCores)
      
      parallel::clusterExport(cl=cl, varlist=c("Sorted", "NewSpan", "NodesPT"),
                              envir = environment())
      parallel::clusterSetRNGStream(cl, iseed = 0L)
      AllFit <- pbapply::pbapply(ExpMat[names(Sorted),], 2, FitFun, cl = cl)
      # AllFit <- parallel::parApply(cl, ExpMat[names(Sorted),], 2, FitFun)
      parallel::stopCluster(cl)
      tictoc::toc()
      
    }
    

    RetVal <- lapply(AllFit, function(x){x$se.fit/x$fit}) %>%
      sapply(., AggFun)

    names(RetVal) <- colnames(ExpMat)

    return(RetVal)

  }

  if(Mode == "LinearOnCircleNodes"){

    print(paste("Computing correlations on", ncol(ExpMat), "genes and", nrow(Edges), "pseudotime segments."))

    NodeOrder <- igraph::subgraph_isomorphisms(target = Net, pattern = igraph::make_lattice(igraph::vcount(Net), circular = TRUE))[[1]] %>%
      names(.)

    Pt <- rpgraph2::getPseudotime(X = ExpMat, Edges = Edges, ProjStruct = ProjStruct, EdgeSeq = c(NodeOrder, NodeOrder[1]))

    Pt <- rpgraph2::getPseudotime(X = ExpMat, Edges = Edges, ProjStruct = ProjStruct, EdgeSeq = c(NodeOrder, NodeOrder[1]))

    ExtPT <- c(Pt$Pt - Pt$PathLen, Pt$Pt, Pt$Pt + Pt$PathLen)

    Selected <- ExtPT >= 0 - Pt$PathLen*.01 & ExtPT <= Pt$PathLen + Pt$PathLen*.01

    X <- ExtPT[Selected]

    AllFit <- apply(ExpMat, 2, function(x){

      Y <- rep(x, 3)[Selected]

      sapply(2:length(Pt$NodePos), function(i){

        ToUse <- (X <= Pt$NodePos[i] & X >= Pt$NodePos[i-1])

        tX <- X[ToUse]
        tY <- Y[ToUse]

        if( sum(ToUse)>3 & length(unique(tX))>1 & length(unique(tY))>1 ){
          cor.test(tY, tX, method = "spe")$p.value
          } else {
          return(1)
        }

      }) %>% AggFun

    })

    names(AllFit) <- colnames(ExpMat)

    return(AllFit)

  }

}

