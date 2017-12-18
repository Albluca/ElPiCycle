#' Fit a circle to the data
#'
#' @param Data
#' @param NonG0Cell
#' @param InitStructNodes
#' @param PlanVarLimitIC
#' @param nNodes
#' @param PlanVarLimit
#'
#' @return
#' @export
#'
#' @examples
FitCircle <- function(Data, NonG0Cell, InitStructNodes, PlanVarLimitIC, PlanVarLimit, nNodes,
                      Phase1_Overwrite = Inf, Phase2_Overwrite = Inf, nCores = 1, PlotIntermediate = FALSE) {

  print("Circle fitting")
  print("Fitting initial circle")

  # Step I - Construct the base circle

  if(!is.null(NonG0Cell)){
    X = Data[rownames(Data) %in% NonG0Cell,]
  } else {
    X = Data
  }

  if(nCores > 1){
    print(paste("Creating a cluster with", nCores, "nodes"))
    cl = parallel::makeCluster(nCores)
    parallel::clusterExport(cl, varlist = c("X"), envir=environment())
  } else {
    cl <- 1
  }

  BasicCircData <- ElPiGraph.R::computeElasticPrincipalCircle(X = X,
                                                           NumNodes = 3, Do_PCA = FALSE,
                                                           CenterData = FALSE, n.cores = cl,
                                                           drawAccuracyComplexity = FALSE,
                                                           drawPCAView = FALSE, drawEnergy = FALSE)


  # BasicCircData[[length(BasicCircData)]][["Phase"]] <- 1

  PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$NodePositions, center = TRUE, scale. = FALSE, retx = TRUE)
  PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

  UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$NodePositions)

  if(is.finite(Phase1_Overwrite)){
    print("Using phase 1 overwrite")
    Cond <- UsedNodes < Phase1_Overwrite
  } else {
    Cond <- UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC
  }

  while(Cond){

    print("Expanding initial circle")

    # Contiune to add node untill the circle remains planar

    BasicCircData <- append(BasicCircData,
                            ElPiGraph.R::computeElasticPrincipalCircle(X = X,
                                                                    NumNodes = UsedNodes + 1,
                                                                    InitNodePositions = BasicCircData[[length(BasicCircData)]]$NodePositions,
                                                                    InitEdges = BasicCircData[[length(BasicCircData)]]$Edges$Edges,
                                                                    Do_PCA = FALSE, CenterData = FALSE,
                                                                    drawAccuracyComplexity = FALSE,
                                                                    n.cores = cl, ComputeMSEP = FALSE,
                                                                    drawPCAView = FALSE, drawEnergy = FALSE)

                            )

    if(PlotIntermediate){
      p <- PlotPG(X = X, TargetPG = BasicCircData[[length(BasicCircData)]])
      Results <- tryCatch(print(p[[1]]), error = function(e){print("Error encontered during plotting. Skipping")})
    }


    # rpgraph2::PlotPG(Data[rownames(Data) %in% NonG0Cell,], BasicCircData[[length(BasicCircData)]], GroupsLab = rep("Non-G0", length(NonG0Cell)))

    BasicCircData[[length(BasicCircData)]][["Phase"]] <- 1

    UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$NodePositions)

    PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$NodePositions, center = TRUE, scale. = FALSE, retx = FALSE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    print(paste("Initial circle: Nodes = ", UsedNodes, "PercPlan = ", PlanPerc))

    if(is.finite(Phase1_Overwrite)){
      print("Using phase 1 overwrite")
      Cond <- UsedNodes < Phase1_Overwrite
    } else {
      Cond <- UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC
    }

  }

  # Step II - Using curves

  if(is.finite(Phase2_Overwrite)){
    print("Using phase 2 overwrite")
    Cond <- UsedNodes < Phase2_Overwrite
  } else {
    Cond <- UsedNodes < nNodes & PlanPerc > PlanVarLimit
  }

  while(Cond){

    print("Extending circle")

    X = Data

    if(nCores > 1){
      parallel::clusterExport(cl, varlist = c("X"), envir=environment())
    }

    BasicCircData <- append(BasicCircData,
                            ElPiGraph.R::computeElasticPrincipalCircle(X = X, NumNodes = UsedNodes + 1,
                                                                    InitNodePositions = BasicCircData[[length(BasicCircData)]]$NodePositions,
                                                                    InitEdges = BasicCircData[[length(BasicCircData)]]$Edges$Edges,
                                                                    Do_PCA = FALSE, CenterData = FALSE,
                                                                    n.cores = cl, ComputeMSEP = FALSE,
                                                                    drawAccuracyComplexity = FALSE,
                                                                    drawPCAView = FALSE, drawEnergy = FALSE)
                            )

    # rpgraph2::PlotPG(Data[rownames(Data) %in% NonG0Cell,], BasicCircData[[length(BasicCircData)]])

    BasicCircData[[length(BasicCircData)]][["Phase"]] <- 2

    # Net <- rpgraph2::ConstructGraph(PrintGraph = BasicCircData[[length(BasicCircData)]])

    PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$NodePositions, center = TRUE, scale. = FALSE, retx = FALSE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)

    UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$NodePositions)

    print(paste("Final circle: Nodes = ", UsedNodes, "PercPlan = ", PlanPerc))

    if(is.finite(Phase2_Overwrite)){
      print("Using phase 2 overwrite")
      Cond <- UsedNodes < Phase2_Overwrite
    } else {
      Cond <- UsedNodes < nNodes & PlanPerc > PlanVarLimit
    }

  }

  if(nCores > 1){
    parallel::stopCluster(cl)
  }

  return(BasicCircData)

}







#'
#'
#'
#' #' Fit a lasso to the data
#' #'
#' #' @param Data
#' #' @param NonG0Cell
#' #' @param InitStructNodes
#' #' @param PlanVarLimitIC
#' #' @param UseTree
#' #' @param ForceLasso
#' #' @param nNodes
#' #' @param MinBranDiff
#' #' @param PlanVarLimit
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' FitLasso <- function(Data, NonG0Cell, InitStructNodes, PlanVarLimitIC, PlanVarLimit, UseTree, ForceLasso, nNodes, MinBranDiff) {
#'
#'
#'   print("Lasso fitting")
#'   print("Fitting initial circle")
#'
#'   # Step I - Construct the base circle
#'
#'   BasicCircData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
#'                                                 Method = 'CircleConfiguration', NodeStep = 1)
#'
#'   PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'   PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'   UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1
#'
#'   while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){
#'
#'     print("Expanding initial circle")
#'
#'     # Contiune to add node untill the circle remains planar
#'
#'     BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
#'                                                                         NumNodes = UsedNodes + 1,
#'                                                                         Method = 'CircleConfiguration',
#'                                                                         NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
#'                                                                         Edges = BasicCircData[[length(BasicCircData)]]$Edges))
#'
#'     UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
#'
#'     PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'     PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'   }
#'
#'   if(UseTree | ForceLasso){
#'
#'     print("Branching initial circle")
#'
#'     # Step II - Construct the initial tree
#'
#'     BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
#'                                                                         NumNodes = UsedNodes + 1,
#'                                                                         Method = 'DefaultPrincipalTreeConfiguration',
#'                                                                         NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
#'                                                                         Edges = BasicCircData[[length(BasicCircData)]]$Edges))
#'
#'
#'     PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'     PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'     UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1
#'
#'     while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
#'
#'       print("Keep Branching")
#'
#'       # Step IIa - keep constructing trees
#'
#'       BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
#'                                                                           NumNodes = UsedNodes + 1,
#'                                                                           Method = 'DefaultPrincipalTreeConfiguration',
#'                                                                           NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
#'                                                                           Edges = BasicCircData[[length(BasicCircData)]]$Edges))
#'
#'
#'       UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
#'
#'       PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'       PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'       Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)
#'       if(max(igraph::degree(Net))>2){
#'         break()
#'       }
#'     }
#'
#'   }
#'
#'   # Step III - Using curves
#'
#'   while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
#'
#'     print("Extending circle and branches")
#'
#'     BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
#'                                                                         NumNodes = nrow(BasicCircData[[length(BasicCircData)]]$Nodes)+1,
#'                                                                         Method = 'CurveConfiguration',
#'                                                                         NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
#'                                                                         Edges = BasicCircData[[length(BasicCircData)]]$Edges))
#'
#'     Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)
#'
#'     PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'     PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'     UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
#'
#'     if(sum(igraph::degree(Net)>2)>1){
#'
#'       print("Multiple branches detected ... trying to select one")
#'
#'       # There are two branhces, this is no good ...
#'       # Look for the biggest circle
#'
#'       CircSize <- 3
#'
#'       while(CircSize < igraph::vcount(Net)){
#'
#'         Template <- igraph::graph.ring(n = CircSize, directed = FALSE, mutual = FALSE, circular = TRUE)
#'
#'         CircSize <- CircSize + 1
#'
#'         CircleInData <- igraph::graph.get.subisomorphisms.vf2(Net, Template)
#'
#'         if(length(CircleInData)>0){
#'           # We found the cicle
#'
#'           EndPoint <- igraph::V(Net)[igraph::degree(Net) == 1]
#'           Branches <- igraph::V(Net)[igraph::degree(Net)>2]
#'           DistMat <- igraph::distances(Net, EndPoint, Branches)
#'
#'           BranchesLen <- apply(DistMat, 1, min)
#'
#'           VertToRemove <- NULL
#'
#'           if(max(BranchesLen[-which.max(BranchesLen)] - max(BranchesLen)) <= -MinBranDiff){
#'
#'             print("Dominating brach found ... pruning the shortest one")
#'
#'             # There is a "dominating"" branch
#'             ToKeep <- names(which.max(BranchesLen))
#'             for(i in 1:nrow(DistMat)){
#'
#'               Source <- rownames(DistMat)[i]
#'               if(ToKeep == Source){
#'                 next()
#'               }
#'
#'               Target <- names(which.min(DistMat[i,]))
#'               VertToRemove <- c(VertToRemove, names(igraph::get.shortest.paths(Net, from = Target, to = Source)$vpath[[1]][-1]))
#'             }
#'
#'             PrunedStruct <- BasicCircData[[length(BasicCircData)]]
#'
#'             NodesToRem <- as.integer(unlist(lapply(strsplit(VertToRemove, "V_"), "[[", 2)))
#'             PrunedStruct$Nodes <- PrunedStruct$Nodes[-NodesToRem, ]
#'             PrunedStruct$Edges <-
#'               PrunedStruct$Edges[!(PrunedStruct$Edges[,1] %in% NodesToRem | PrunedStruct$Edges[,2] %in% NodesToRem), ]
#'
#'             for(VerVal in 1:max(PrunedStruct$Edges)){
#'               if(any(PrunedStruct$Edges == VerVal)){
#'                 next()
#'               } else {
#'                 PrunedStruct$Edges[PrunedStruct$Edges>VerVal] <-
#'                   PrunedStruct$Edges[PrunedStruct$Edges>VerVal] - 1
#'               }
#'
#'             }
#'
#'             PrunedStruct$Method <- "ManualPruning"
#'
#'             BasicCircData[[length(BasicCircData) + 1]] <- PrunedStruct
#'
#'             PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'             PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'             UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
#'
#'             break()
#'
#'           }
#'
#'
#'         }
#'
#'
#'       }
#'
#'     }
#'
#'   }
#'
#'   return(BasicCircData)
#'
#' }
#'
#'
#'
#'
#'
#'
#' #' Fit a line to the data
#' #'
#' #' @param Data
#' #' @param NonG0Cell
#' #' @param InitStructNodes
#' #' @param PlanVarLimitIC
#' #' @param nNodes
#' #' @param PlanVarLimit
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' FitLine <- function(Data, NonG0Cell, InitStructNodes, PlanVarLimitIC, PlanVarLimit, nNodes) {
#'
#'   print("Line fitting")
#'   print("Fitting initial line")
#'
#'   # Step I - Construct the base circle
#'
#'   BasicLineData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
#'                                                 Method = 'CurveConfiguration', NodeStep = 1)
#'
#'   PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'   PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'   UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes) - 1
#'
#'   while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){
#'
#'     print("Expanding initial line")
#'
#'     # Contiune to add node untill the circle remains planar
#'
#'     BasicLineData <- append(BasicLineData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
#'                                                                         NumNodes = UsedNodes + 1,
#'                                                                         Method = 'CircleConfiguration',
#'                                                                         NodesPositions = BasicLineData[[length(BasicLineData)]]$Nodes,
#'                                                                         Edges = BasicLineData[[length(BasicLineData)]]$Edges))
#'
#'
#'     UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes)
#'
#'     PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'     PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'   }
#'
#'   # Step II - Using curves
#'
#'   while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
#'
#'     print("Extending line")
#'
#'     BasicLineData <- append(BasicLineData, computeElasticPrincipalGraph(Data = Data,
#'                                                                         NumNodes = nrow(BasicLineData[[length(BasicLineData)]]$Nodes)+1,
#'                                                                         Method = 'CurveConfiguration',
#'                                                                         NodesPositions = BasicLineData[[length(BasicLineData)]]$Nodes,
#'                                                                         Edges = BasicLineData[[length(BasicLineData)]]$Edges))
#'
#'
#'     Net <- ConstructGraph(Results = BasicLineData[[length(BasicLineData)]], DirectionMat = NULL, Thr = NULL)
#'
#'     PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
#'     PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
#'
#'     UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes)
#'
#'   }
#'
#'   return(BasicLineData)
#'
#' }
#'
#'
#'
#'
