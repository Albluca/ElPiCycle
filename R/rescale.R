#' Title
#'
#' @param OrderedData
#' @param RescaleVal
#'
#' @return
#' @export
#'
#' @examples
RescaleStagePT <- function(OrderedData,
                           RescaleVal = list("G0" = c(0, .25),
                                             "G1" = c(.25, .5),
                                             "S" = c(.5, .75),
                                             "G2" = c(.75, 1))) {

  Ranges <- lapply(OrderedData$StageOnNodes, range)

  CatNodes <- lapply(1:length(RescaleVal), function(i){
    names(RescaleVal) %>%
      "[[" (i) %>%
      grep(., names(OrderedData$StageOnNodes), fixed = TRUE)
  })

  names(CatNodes) <- names(RescaleVal)

  CatNodes[[1]] <- setdiff(CatNodes[[1]], CatNodes[[length(CatNodes)]])

  RescaledPT <- 0

  for(i in 1:(length(CatNodes)-1)){
    OverLapCat <- intersect(CatNodes[[i]], CatNodes[[i+1]])

    NoOverlap <- FALSE

    if(length(OverLapCat)>0){
      if(
        lapply(OrderedData$StageOnNodes[OverLapCat], length) %>%
        unlist %>%
        ">" (0)
      ){
        OverLapnodes <- OrderedData$StageOnNodes[OverLapCat]
        ChangePoint <- OrderedData$NodesPT[unlist(OverLapnodes)] %>%
          range %>%
          mean
      } else {
        NoOverlap <- TRUE
      }
    }

    if(NoOverlap | length(OverLapCat) == 0){

      LowPoint <- CatNodes[[i]] %>% unlist(.) %>%
        OrderedData$StageOnNodes[.] %>% unlist(.) %>% max(.)

      HighPoint <- CatNodes[[i+1]] %>% unlist(.) %>%
        OrderedData$StageOnNodes[.] %>% unlist(.) %>% min(.)

      if(is.infinite(HighPoint) & !is.infinite(LowPoint)){
        ChangePoint <- OrderedData$NodesPT[LowPoint] %>% range %>% mean
      }

      if(!is.infinite(HighPoint) & is.infinite(LowPoint)){
        ChangePoint <- OrderedData$NodesPT[HighPoint] %>% range %>% mean
      }

      if(!is.infinite(HighPoint) & !is.infinite(LowPoint)){
        ChangePoint <- OrderedData$NodesPT[c(LowPoint, HighPoint)] %>% range %>% mean
      }

      if(is.infinite(HighPoint) & is.infinite(LowPoint)){
        ChangePoint <- RescaledPT[length(RescaledPT)]
      }

    }

    RescaledPT <- c(RescaledPT, ChangePoint)

  }

  RescaledPT <- c(RescaledPT, max(OrderedData$NodesPT))

  MatricixedList <- sapply(RescaleVal, rbind)

  tFun <- approxfun(
    x = RescaledPT,
    y = c(MatricixedList[1,], MatricixedList[2,ncol(MatricixedList)]),
    method = "linear"
    )

  tVect <- tFun(OrderedData$CellsPT)
  names(tVect) <- names(OrderedData$CellsPT)
  OrderedData$CellsPT <- tVect

  tVect <- tFun(OrderedData$NodesPT)
  names(tVect) <- names(OrderedData$NodesPT)
  OrderedData$NodesPT <- tVect

  tVect <- t(apply(OrderedData$RecCoord[,1:3], 1, tFun))
  OrderedData$RecCoord[,1:3] <- tVect

  return(OrderedData)

}
