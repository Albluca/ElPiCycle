
################################################################################
#
# Function used to project cells on a circle or lasso ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param nNodes
#' @param GraphType
#' @param PlanVarLimit
#' @param PlanVarLimitIC
#' @param InitStructNodes
#' @param DataMat 
#' @param nCores 
#' @param NonG0Cell
#'
#' @return
#' @export
#'
#' @examples
#'
ProjectAndCompute <- function(DataMat,
                              nNodes = 30,
                              GraphType = 'Circle',
                              PlanVarLimit = .9,
                              PlanVarLimitIC = NULL,
                              InitStructNodes = 15,
                              nCores = 1,
                              NonG0Cell){

  # if PlanVarLimitIC is not set it will be the same as PlanVarLimit
  if(is.null(PlanVarLimitIC)){
    PlanVarLimitIC <- PlanVarLimit
  }

  # Fitting a circle
  if(GraphType == 'Circle') {
    FitData <- FitCircle(Data = DataMat, NonG0Cell =  NonG0Cell, InitStructNodes =  InitStructNodes,
                         PlanVarLimitIC = PlanVarLimitIC, PlanVarLimit = PlanVarLimit, nNodes = nNodes, nCores = nCores)
  }

  return(list(FitData = FitData, GraphType = GraphType, PlanVarLimit = PlanVarLimit,
              PlanVarLimitIC = PlanVarLimitIC, InitStructNodes = InitStructNodes, nNodes = nNodes))
}












################################################################################
#
# Function used to construct and optimize the cycle ------------------------------------------------------
#
################################################################################



#' Construct and optimize cycle
#'
#' @param DataSet
#' @param StartSet
#' @param VarThr
#' @param nNodes
#' @param Categories
#' @param GraphType
#' @param PlanVarLimit
#' @param PlanVarLimitIC
#' @param InitStructNodes
#' @param EstProlif
#' @param QuaThr
#' @param NonG0Cell
#' @param MinProlCells
#' @param AddGenePerc
#' @param SelThr1
#' @param SelThr2
#' @param MadsThr
#' @param OutThr.Gene_Expression
#' @param OutThr.Gene_Count
#' @param OutThr.Gene_Space
#' @param LogSpace
#' @param Do.LogPC
#' @param Do.PCA
#' @param Center.PCA
#' @param GeneSelMode
#' @param SelGeneThr
#' @param SelGeneAggFun
#' @param Span
#' @param nCores
#' @param Title
#'
#' @return
#' @export
#'
#' @examples
SelectGenesOnGraph <- function(
  # Base parameters
  DataSet,
  StartSet,
  Categories = NULL,
  # Expression filtering
  OutThr.Gene_Expression = 3,
  OutThr.Gene_Count = 3,
  OutThr.Gene_Space = 3,
  LogSpace = TRUE,
  # Proliferation estimation
  NonG0Cell = NULL,
  EstProlif = "MeanPerc",
  QuaThr = .5,
  MinProlCells = 50,
  # Transformation (Log, PCA)
  Do.LogPC = TRUE,
  Do.PCA = TRUE,
  Center.PCA = FALSE,
  VarThr = .99,
  # Principal circle estimation
  nNodes = 40,
  GraphType = "Circle",
  PlanVarLimit = .85,
  PlanVarLimitIC = .9,
  InitStructNodes = 20,
  # Gene selection
  GeneSelMode = "SmoothOnCircleNodes",
  AddGenePerc = 5,
  SelThr1 = .95,
  SelThr2 = .99,
  MadsThr =  1,
  SelGeneThr = NULL,
  SelGeneAggFun = median,
  Span = .5,
  # Parallelization
  nCores = 1,
  # Plot labels
  Title = ''
  ) {

  # Initialize the structures
  Steps <- list()
  UsedGenes <- list()

  # Filter data - Remove non detected genes
  print(paste(sum(colSums(DataSet)>0), "genes expressed in", nrow(DataSet), "cells"))
  SelGenes <- colSums(DataSet)>0

  # Keep track of the genes being used
  UsedGenes[[1]] <- intersect(StartSet, colnames(DataSet[,SelGenes]))

  # inizialize Categories if not specified
  if(is.null(Categories)){
    Categories <- rep("N/A", nrow(DataSet))
  }
  
  # Make Categories a factor if it is not not specified
  if(!is.factor(Categories)){
    print("Categories will be transformed in a factor")
    Categories <- factor(Categories)
  }

  if(length(UsedGenes[[1]]) < 10){
    stop("Number of selected genes < 10. Impossible to proceed")
  }

  # Filter data - Remove outlier cells
  ToFilter <- FilterData(DataMat = DataSet[, SelGenes],
                         OutThr.Gene_Expression = OutThr.Gene_Expression,
                         OutThr.Gene_Count = OutThr.Gene_Count,
                         OutThr.Gene_Space = OutThr.Gene_Space,
                         LogSpace = LogSpace)
  print(paste(sum(!ToFilter), "cells passed filtering"))

  # Subset the matrix and categories to consider only the genes used
  FiltExpmat <- DataSet[!ToFilter, UsedGenes[[1]]]
  FiltCat <- Categories[!ToFilter]
  print(paste(length(UsedGenes[[1]]), "genes selected in", nrow(DataSet), "cells"))

  # Transform data
  if(Do.LogPC){
    print("Using pseudocount (log10(x+1))")
    FiltExpmat <- log10(FiltExpmat + 1)
  }

  # Filter data - Select proliferative cells
  if(is.null(NonG0Cell)){
    ProlCells <- EstimateProliferativeCell(DataMat = FiltExpmat, EstProlif = EstProlif, QuaThr = QuaThr)
    NonG0Cell <- ProlCells$NonG0Cell
  }

  if(length(NonG0Cell) < MinProlCells){
    print("Not enought proliferative cells detected. Ignoring them")
    NonG0Cell <- NULL
  }

  if(Do.PCA){
    PCAData <- prcomp(FiltExpmat, retx = TRUE, scale. = FALSE, center = Center.PCA)
    FiltExpmat <- PCAData$x
  }

  Steps[[1]] <- ProjectAndCompute(DataMat = FiltExpmat,
                                  nNodes = nNodes,
                                  GraphType = GraphType,
                                  PlanVarLimit = PlanVarLimit,
                                  PlanVarLimitIC = PlanVarLimitIC,
                                  InitStructNodes = InitStructNodes,
                                  NonG0Cell = NonG0Cell, nCores = nCores)


  p <- ElPiGraph.R::PlotPG(FiltExpmat, TargetPG = Steps[[1]]$FitData[[length(Steps[[1]]$FitData)]], GroupsLab = FiltCat, PGCol = "black",
                   Main = paste0(Title, " / Step 1 / ", length(UsedGenes[[1]]), " genes "))

  Results <- tryCatch(print(p[[1]]), error = function(e){print("Error encontered during plotting. Skipping")})


  # Select the genes that are closer to the original curve

  CONVERGED <- FALSE
  i = 2

  while(!CONVERGED){

    print("Phase I - Contraction")

    Partition <- ElPiGraph.R::PartitionData(X = FiltExpmat,
                                         NodePositions = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$NodePositions,
                                         SquaredX = rowSums(FiltExpmat^2)
                                         )
    Net <- ElPiGraph.R::ConstructGraph(PrintGraph = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]])

    ProjStruct <- ElPiGraph.R::project_point_onto_graph(X = FiltExpmat,
                                       NodePositions = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$NodePositions,
                                       Edges = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$Edges$Edges,
                                       Partition = Partition$Partition)

    print("Selecting genes for the restriction set")

    tictoc::tic()
    if(Do.LogPC){
      GenesThr <- SelectGenes(Partition = Partition$Partition,
                              Net = Net,
                              ExpMat = log10(DataSet[!ToFilter, UsedGenes[[i-1]]] + 1),
                              Mode = GeneSelMode,
                              AggFun = SelGeneAggFun,
                              Span = Span,
                              Edges = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$Edges$Edges,
                              ProjStruct = ProjStruct,
                              nCores = nCores)
    } else {
      GenesThr <- SelectGenes(Partition = Partition$Partition,
                              Net = Net,
                              ExpMat = DataSet[!ToFilter, UsedGenes[[i-1]]],
                              Mode = GeneSelMode,
                              AggFun = SelGeneAggFun,
                              Span = Span,
                              Edges = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$Edges$Edges,
                              ProjStruct = ProjStruct,
                              nCores = nCores)
    }
    tictoc::toc()

    GenesThr <- GenesThr[is.finite(GenesThr)]

    boxplot(GenesThr)

    if(i == 2 & is.null(SelGeneThr)){
      SelGeneThr <- median(GenesThr) + MadsThr*mad(GenesThr)
    }

    abline(h= SelGeneThr)

    UsedGenes[[i]] <- names(which(GenesThr < SelGeneThr))

    print(
      paste(
        length(setdiff(UsedGenes[[length(UsedGenes) - 1]], UsedGenes[[length(UsedGenes)]])),
        "genes removed /",
        length(UsedGenes[[length(UsedGenes)]]),
        "genes are now selected"
        )
      )

    if(length(UsedGenes[[length(UsedGenes)]]) == 0){
      print("No gene selected. Impossible to proceed!")
      return(NULL)
    }

    # Subset the matrix to consider only the genes used
    FiltExpmat <- DataSet[!ToFilter, UsedGenes[[i]]]
    print(paste(length(UsedGenes[[i]]), "genes selected in", nrow(DataSet), "cells"))

    # Filter data - Select proliferative cells
    if(is.null(NonG0Cell)){
      ProlCells <- EstimateProliferativeCell(DataMat = FiltExpmat, EstProlif = EstProlif, QuaThr = QuaThr)
      NonG0Cell <- ProlCells$NonG0Cell
    }

    if(length(NonG0Cell) < MinProlCells){
      print("Not enought proliferative cells detected. Ignoring them")
      NonG0Cell <- NULL
    }

    # Transform data

    if(Do.LogPC){
      FiltExpmat <- log10(FiltExpmat + 1)
    }

    if(Do.PCA){
      PCAData <- prcomp(FiltExpmat, retx = TRUE, scale. = FALSE, center = Center.PCA)
      FiltExpmat <- PCAData$x
    }

    Steps[[i]] <- ProjectAndCompute(DataMat = FiltExpmat,
                                   nNodes = nNodes,
                                   GraphType = GraphType,
                                   PlanVarLimit = PlanVarLimit,
                                   PlanVarLimitIC = PlanVarLimitIC,
                                   InitStructNodes = InitStructNodes,
                                   NonG0Cell = NonG0Cell,
                                   nCores = nCores)

    p <- ElPiGraph.R::PlotPG(FiltExpmat, TargetPG = Steps[[i]]$FitData[[length(Steps[[i]]$FitData)]], GroupsLab = FiltCat,
                     PGCol = "EpG",
                     Main = paste0(Title, " / Step ", i, " / ", length(UsedGenes[[i]]), " genes "))

    Results <- tryCatch(print(p[[1]]), error = function(e){print("Error encontered during plotting. Skipping")})

    i = i + 1

    if(length(UsedGenes[[i-1]])/length(UsedGenes[[i - 2]]) > SelThr1){
      CONVERGED <- TRUE
    }

  }

  StageIEnd <- i

  # Extend the geneset by including other genes that are closer to the original circle

  CONVERGED2 <- FALSE

  if(AddGenePerc == 0){
    CONVERGED2 <- TRUE
  }
  
  while(!CONVERGED2){

    ######### WORK HERE

    print("Phase II - Expansion")

    Partition <- ElPiGraph.R::PartitionData(X = FiltExpmat,
                                         NodePositions = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$NodePositions,
                                         SquaredX = rowSums(FiltExpmat^2)
    )
    Net <- ElPiGraph.R::ConstructGraph(PrintGraph = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]])

    ProjStruct <- ElPiGraph.R::project_point_onto_graph(X = FiltExpmat,
                                                     NodePositions = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$NodePositions,
                                                     Edges = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$Edges$Edges,
                                                     Partition = Partition$Partition)

    print("Selecting genes for the expansion set")

    tictoc::tic()
    if(Do.LogPC){
      GenesThr <- SelectGenes(Partition = Partition$Partition,
                              Net = Net,
                              ExpMat = log10(DataSet[!ToFilter, ] + 1),
                              Mode = GeneSelMode,
                              AggFun = SelGeneAggFun,
                              Span = Span,
                              Edges = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$Edges$Edges,
                              ProjStruct = ProjStruct,
                              nCores = nCores)
    } else {
      GenesThr <- SelectGenes(Partition = Partition$Partition,
                              Net = Net,
                              ExpMat = DataSet[!ToFilter, ],
                              Mode = GeneSelMode,
                              AggFun = SelGeneAggFun,
                              Span = Span,
                              Edges = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]]$Edges$Edges,
                              ProjStruct = ProjStruct,
                              nCores = nCores)
    }
    tictoc::toc()

    GenesThr <- GenesThr[is.finite(GenesThr)]

    AddGenes <- sort(GenesThr, decreasing = FALSE)[1:round(AddGenePerc*length(UsedGenes[[i-1]])/100)]
    AddGenes <- AddGenes[AddGenes < SelGeneThr]
    UsedGenes[[i]] <- union(names(AddGenes), UsedGenes[[i-1]])


    print(
      paste(
        length(setdiff(UsedGenes[[length(UsedGenes)]], UsedGenes[[length(UsedGenes) - 1]])),
        "genes added /",
        length(UsedGenes[[length(UsedGenes)]]),
        "genes are now selected"
      )
    )

    if(length(UsedGenes[[length(UsedGenes)]]) == 0){
      print("No gene selected. Impossible to proceed!")
      return(NULL)
    }

    # Subset the matrix to consider only the genes used
    FiltExpmat <- DataSet[!ToFilter, UsedGenes[[i]]]
    print(paste(length(UsedGenes[[i]]), "genes selected in", nrow(DataSet), "cells"))

    # Filter data - Select proliferative cells
    if(is.null(NonG0Cell)){
      ProlCells <- EstimateProliferativeCell(DataMat = FiltExpmat, EstProlif = EstProlif, QuaThr = QuaThr)
      NonG0Cell <- ProlCells$NonG0Cell
    }

    if(length(NonG0Cell) < MinProlCells){
      print("Not enought proliferative cells detected. Ignoring them")
      NonG0Cell <- NULL
    }

    # Transform data if needed

    if(Do.LogPC){
      FiltExpmat <- log10(FiltExpmat + 1)
    }

    if(Do.PCA){
      PCAData <- prcomp(FiltExpmat, retx = TRUE, scale. = FALSE, center = Center.PCA)
      FiltExpmat <- PCAData$x
    }

    Steps[[i]] <- ProjectAndCompute(DataMat = FiltExpmat,
                                    nNodes = nNodes,
                                    GraphType = GraphType,
                                    PlanVarLimit = PlanVarLimit,
                                    PlanVarLimitIC = PlanVarLimitIC,
                                    InitStructNodes = InitStructNodes,
                                    NonG0Cell = NonG0Cell,
                                    nCores = nCores)

    p <- ElPiGraph.R::PlotPG(FiltExpmat, TargetPG = Steps[[i]]$FitData[[length(Steps[[i]]$FitData)]], GroupsLab = FiltCat,
                          PGCol = "EpG",
                          Main = paste0(Title, " / Step ", i, " / ", length(UsedGenes[[i]]), " genes "))

    Results <- tryCatch(print(p[[1]]), error = function(e){print("Error encontered during plotting. Skipping")})

    i = i + 1

    if(length(UsedGenes[[i-2]])/length(UsedGenes[[i-1]]) > SelThr2){
      CONVERGED2 <- TRUE
    }

  }

  plot(unlist(lapply(UsedGenes, length)))
  abline(v=StageIEnd-.5)
  abline(v=1.5)

  if(Do.PCA){
    return(list(Genes = UsedGenes, PGStructs = Steps,
                FinalStruct = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]],
                FinalExpMat = FiltExpmat, FinalGroup = FiltCat,
                PCARotation = PCAData$rotation, PCACenter = PCAData$center))
  } else {
    return(list(Genes = UsedGenes, PGStructs = Steps,
                FinalStruct = Steps[[i-1]]$FitData[[length(Steps[[i-1]]$FitData)]],
                FinalExpMat = FiltExpmat, FinalGroup = FiltCat,
                PCARotation = NULL, PCACenter = NULL))
  }


}








