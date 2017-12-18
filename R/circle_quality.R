#' Assess quality of the fitting via correlation
#'
#' @param TargetStruct
#' @param ProcStruct
#' @param DistVect
#' @param CorMet
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
DistanceGraphCor <- function(TargetStruct, ProcStruct, DistVect = 4:10, CorMet = "spe") {

  AllDist <- distances(TargetStruct$Net[[length(TargetStruct$Net)]])

  AllDist_Reord <- colnames(AllDist) %>%
    strsplit(., split = '_', fixed = TRUE) %>%
    sapply(., "[[", 2) %>%
    as.integer %>%
    order %>%
    AllDist[., .]

  AllPos <- rbind(TargetStruct$Data, TargetStruct$PrinGraph$Nodes) %>%
    stats::dist(.) %>%
    as.matrix(.)

  CellDist <- AllPos[1:nrow(TargetStruct$Data), -c(1:nrow(TargetStruct$Data))]

  AssignedNodes <- apply(CellDist, 1, which.min)

  CorVal <- sapply(DistVect, function(thr){
    sapply(X = 1:length(AssignedNodes), FUN = function(i){
      SelDistVect <- AllDist_Reord[AssignedNodes[i],]
      cor(SelDistVect[SelDistVect <= thr], CellDist[i, SelDistVect <= thr], method = CorMet)
    })
  })

  rownames(CorVal) <- names(AssignedNodes)
  colnames(CorVal) <- paste("P", DistVect, sep = "=")

  return(CorVal)

}
