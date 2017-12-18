#' Title
#'
#' @param TargetStruct
#' @param GraphType
#' @param InitStructNodes
#' @param PlanVarLimit
#' @param PlanVarLimitIC
#' @param nNodes
#' @param ForceLasso
#' @param MinBranDiff
#' @param nSamples
#' @param CellPerc
#' @param MinProlCells
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
BootStrapCurve <- function(
  TargetStruct,
  InitStructNodes = 15,
  PlanVarLimit = .85,
  PlanVarLimitIC = .9,
  nNodes = 40,
  ForceLasso = FALSE,
  MinBranDiff = 2,
  nSamples = 10,
  CellPerc = .95,
  MinProlCells = 10,
  VarThr = 1,
  Cores = 1) {

  GraphType = 'Circle'

  # TargetStruct <- Trap_D0.Data$Analysis$PGStructs[[length(Trap_D0.Data$Analysis$PGStructs)]]
  ToSample <- round(CellPerc*nrow(TargetStruct$Data))

  if(is.null(PlanVarLimitIC)){
    PlanVarLimitIC <- PlanVarLimit + .5*(1-PlanVarLimit)
  }

  SampledStruct <- list()

  PhaseMax <- lapply(TargetStruct$IntGrahs, "[[", "Phase") %>%
    unlist %>%
    aggregate(x=1:length(.), by=list(.), FUN=max)

  NodesByGraph <- lapply(TargetStruct$IntGrahs, "[[", "Nodes") %>%
    lapply(., nrow) %>%
    unlist

  Phase1_Nodes <- NodesByGraph[PhaseMax[1, 2]]
  Phase2_Nodes <- NodesByGraph[PhaseMax[2, 2]]


  for(i in 1:nSamples){

    print(paste(rep("#", 20+4), collapse = ""))
    print(paste("Sample", i))
    print(paste(rep("#", 20+4), collapse = ""))

    nCol <- (cumsum(TargetStruct$ExpVar) <= VarThr) %>%
      which %>%
      max %>%
      min(., ncol(TargetStruct$Data))

    ExpData <- TargetStruct$Data[sample(1:nrow(TargetStruct$Data), ToSample), 1:nCol]
    SelNonG0Cell <- intersect(TargetStruct$NonG0Cell,
                              rownames(ExpData))




    if(length(SelNonG0Cell)>MinProlCells){

      print("Using strongly proliferative cells for initial fitting")
      UseTree <- TRUE

    } else {

      print("Unable to find a sufficent number of strongly proliferative cells")
      SelNonG0Cell <- rownames(ExpData)
      UseTree <- FALSE

    }


    if(GraphType == 'Lasso') {
      FitData <- FitLasso(Data = ExpData,
                          NonG0Cell = SelNonG0Cell,
                          InitStructNodes = InitStructNodes,
                          PlanVarLimitIC = PlanVarLimitIC,
                          PlanVarLimit = PlanVarLimit,
                          UseTree = UseTree,
                          ForceLasso = ForceLasso,
                          nNodes = nNodes,
                          MinBranDiff = MinBranDiff)
    }


    if(GraphType == 'Circle') {
      FitData <- FitCircle(Data = ExpData,
                           NonG0Cell = SelNonG0Cell,
                           InitStructNodes = InitStructNodes,
                           PlanVarLimitIC = PlanVarLimitIC,
                           PlanVarLimit = PlanVarLimit,
                           nNodes = nNodes,
                           Phase1_Overwrite = Phase1_Nodes,
                           Phase2_Overwrite = Phase2_Nodes)
    }


    if(GraphType == 'Line') {
      FitData <- FitLine(Data = ExpData,
                         NonG0Cell = SelNonG0Cell,
                         InitStructNodes = InitStructNodes,
                         PlanVarLimitIC = PlanVarLimitIC,
                         PlanVarLimit = PlanVarLimit,
                         nNodes = nNodes)
    }

    SampledStruct[[i]] <- FitData[[length(FitData)]]

  }


  return(SampledStruct)

}






#' Plot bootstrapped data
#'
#' @param SampledStruct
#' @param PoitMult
#' @param TargetStruct
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
PlotBoot <- function(SampledStruct, TargetStruct, PoitMult = 1.5, Title = "") {

  # PCA of the reference principal graphs

  BasePCA <- prcomp(TargetStruct$PrinGraph$Nodes, retx = TRUE, center = TRUE, scale. = FALSE)

  RefPoints12 <- data.frame(BasePCA$x[,1:2])
  RefPoints13 <- data.frame(BasePCA$x[,c(1,3)])

  RefEdges <- TargetStruct$PrinGraph$Edges


  # Get all of the nodes of the bootstrapped structures and perform (Circle only atm)

  AllBootNodes <- lapply(SampledStruct, "[[", "Nodes") %>%
    do.call(rbind, .)


  MeanBoot <- computeElasticPrincipalGraph(Data = AllBootNodes,
                               NumNodes = PoitMult*(TargetStruct$PrinGraph$Nodes %>% nrow),
                               Method = 'CircleConfiguration')


  # Rotate the nodes of the mean bootstrapped and bootstrapped structures using the PCA of the refernce curve

  RotateMeanBootNodes <- t(t(MeanBoot[[1]]$Nodes) - BasePCA$center[1:ncol(MeanBoot[[1]]$Nodes)]) %>%
    '%*%'(BasePCA$rotation[1:ncol(MeanBoot[[1]]$Nodes), ])

  MeanBootPoint_12 <- data.frame(RotateMeanBootNodes[,1:2], Sample = "mean")
  MeanBootPoint_13 <- data.frame(RotateMeanBootNodes[,c(1,3)], Sample = "mean")

  # Create accessory structures for the structures

  ListRotated <- lapply(SampledStruct, "[[", "Nodes") %>%
    lapply(., function(x) {t(t(x) - BasePCA$center[1:ncol(MeanBoot[[1]]$Nodes)])}) %>%
    lapply(., function(x) {x %*% BasePCA$rotation[1:ncol(MeanBoot[[1]]$Nodes),]})

  PointList_12 <- lapply(1:length(ListRotated), function(i){data.frame(ListRotated[[i]][,1:2], Sample = i)})
  PointList_13 <- lapply(1:length(ListRotated), function(i){data.frame(ListRotated[[i]][,c(1,3)], Sample = i)})

  EdgesList_12 <- lapply(1:length(ListRotated), function(i){
    apply(SampledStruct[[i]]$Edges, 1, function(x){
      data.frame(PC1_1 = ListRotated[[i]][x[1],1],
                 PC2_1 = ListRotated[[i]][x[1],2],
                 PC1_2 = ListRotated[[i]][x[2],1],
                 PC2_2 = ListRotated[[i]][x[2],2],
                 Sample = i)
    })
  }) %>% lapply(., do.call, what = rbind) %>%
    do.call(rbind, .)


  EdgesList_13 <- lapply(1:length(ListRotated), function(i){
    apply(SampledStruct[[i]]$Edges, 1, function(x){
      data.frame(PC1_1 = ListRotated[[i]][x[1],1],
                 PC3_1 = ListRotated[[i]][x[1],3],
                 PC1_2 = ListRotated[[i]][x[2],1],
                 PC3_2 = ListRotated[[i]][x[2],3],
                 Sample = i)
    })
  }) %>% lapply(., do.call, what = rbind) %>%
    do.call(rbind, .)




  SampledPoints <- nrow(PointList_12[[1]])
  SampledEdges <- nrow(EdgesList_12[[1]])




  # Create the base plot using the cellsas reference
  p1 <- ggplot2::ggplot(data = data.frame(t(t(TargetStruct$Data) - BasePCA$center) %*% BasePCA$rotation[,1:2],
                                          Cat = TargetStruct$Categories),
                       mapping = ggplot2::aes(x = PC1, y = PC2, color = Cat)) +

    # Add the points of the bootstrapped principal curve
    # lapply(1:length(PointList_12), function(i){
    #   ggplot2::geom_point(data = data.frame(PointList_12[[i]]),
    #                       mapping = ggplot2::aes(x = PC1, y = PC2, color = as.character(Sample)), inherit.aes = FALSE)
    # }) +

    # Add the cells
    ggplot2::geom_point() +

    # Add the edges of the bootstrapped principal curves
    ggplot2::geom_segment(data = EdgesList_12, ggplot2::aes(x=PC1_1, xend=PC1_2,
                                                            y=PC2_1, yend=PC2_2),
                          alpha = .2, color = "black") +

    # Add the edges of the mean bootstrapped principal curves
    apply(MeanBoot[[1]]$Edges, 1, function(x){data.frame(MeanBootPoint_12[x,])}) %>%
    lapply(., function(x){
      data.frame(PC1_1 = x[1,1], PC1_2 = x[2,1],
                 PC2_1 = x[1,2], PC2_2 = x[2,2], Sample = "mean")
    }) %>%
    do.call(rbind, .) %>%
    ggplot2::geom_segment(data = ., mapping = ggplot2::aes(x = PC1_1, xend = PC1_2,
                                                           y = PC2_1, yend = PC2_2), inherit.aes = FALSE,
                          color = "red", linetype = 'dashed', size = 1.5) +

    # Add the principal curve
    ggplot2::geom_point(data = RefPoints12,
                        mapping = ggplot2::aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "blue") +

    # Add the edges of the principal curve
    apply(RefEdges, 1, function(x){data.frame(RefPoints12[x,])}) %>%
    lapply(., function(x){
      data.frame(PC1_1 = x[1,1], PC1_2 = x[2,1],
                 PC2_1 = x[1,2], PC2_2 = x[2,2])
    }) %>%
    do.call(rbind, .) %>%
    ggplot2::geom_segment(data = ., mapping = ggplot2::aes(x = PC1_1, xend = PC1_2,
                                                           y = PC2_1, yend = PC2_2), inherit.aes = FALSE,
                          color = "blue", size = 1) +

    # Set some aesthetic
    ggplot2::labs(title = Title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


  print(p1)







  # Create the base plot with the cells
  p2 <- ggplot2::ggplot(data = data.frame(t(t(TargetStruct$Data) - BasePCA$center) %*% BasePCA$rotation[,c(1,3)],
                                          Cat = TargetStruct$Categories),
                       mapping = ggplot2::aes(x = PC1, y = PC3, color = Cat)) +

    # Add the points of the bootstrapped principal curve
    # lapply(1:length(PointList_23), function(i){
    #   ggplot2::geom_point(data = data.frame(PointList_23[[i]]),
    #                       mapping = ggplot2::aes(x = PC2, y = PC3, color = as.character(Sample)), inherit.aes = FALSE)
    # }) +

    # Add the cells
    ggplot2::geom_point() +

    # Add the edges of the bootstrapped principal curve
    ggplot2::geom_segment(data = EdgesList_13, ggplot2::aes(x=PC1_1, xend=PC1_2,
                                                            y=PC3_1, yend=PC3_2),
                          alpha = .2, color = "black") +

    # Add the edges of the mean bootstrapped principal curves
    apply(MeanBoot[[1]]$Edges, 1, function(x){data.frame(MeanBootPoint_13[x,])}) %>%
    lapply(., function(x){
      data.frame(PC1_1 = x[1,1], PC1_2 = x[2,1],
                 PC3_1 = x[1,2], PC3_2 = x[2,2], Sample = "mean")
    }) %>%
    do.call(rbind, .) %>%
    ggplot2::geom_segment(data = ., mapping = ggplot2::aes(x = PC1_1, xend = PC1_2,
                                                           y = PC3_1, yend = PC3_2), inherit.aes = FALSE,
                          color = "red", linetype = 'dashed', size = 1.5) +

    # Add the principal curve
    ggplot2::geom_point(data = data.frame(BasePCA$x[,c(1,3)]),
                        mapping = ggplot2::aes(x = PC1, y = PC3), inherit.aes = FALSE, color = "blue") +

    # Add the edges of the principal curve
    apply(RefEdges, 1, function(x){data.frame(RefPoints13[x,])}) %>%
    lapply(., function(x){
      data.frame(PC1_1 = x[1,1], PC1_2 = x[2,1],
                 PC3_1 = x[1,2], PC3_2 = x[2,2])
    }) %>%
    do.call(rbind, .) %>%
    ggplot2::geom_segment(data = ., mapping = ggplot2::aes(x = PC1_1, xend = PC1_2,
                                                           y = PC3_1, yend = PC3_2), inherit.aes = FALSE,
                          color = "blue", size = 1) +

    # Set some aesthetic
    ggplot2::labs(title = Title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


  print(p2)


  return(list("PC12" = p1, "PC13" = p2, MeanBoot = MeanBoot))

}
