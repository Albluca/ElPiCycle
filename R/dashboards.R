#' Title
#'
#' @param InputList
#'
#' @return
#' @export
#'
#' @examples
CompareAcrossData <- function(InputList, CatOrder = NULL) {

  library(shiny)

  nSpans = length(InputList)

  for(i in 1:length(InputList)){

    WorkStruct <- InputList[[i]]$OrderedData

    ValidCells <- intersect(rownames(InputList[[i]]$Expression),
              intersect(rownames(WorkStruct$CellExp),
                        names(WorkStruct$CellsPT))
              )

    InputList[[i]]$Expression <- InputList[[i]]$Expression[rownames(InputList[[i]]$Expression) %in% ValidCells, ]
    WorkStruct$CellExp <- WorkStruct$CellExp[rownames(WorkStruct$CellExp) %in% ValidCells, ]
    WorkStruct$CellsPT <- WorkStruct$CellsPT[names(WorkStruct$CellsPT) %in% ValidCells]

    InputList[[i]]$OrderedData <- WorkStruct

  }

  PlotH <- paste0(length(InputList)*400, "px")


  # Define UI
  ui <- fluidPage(

    # Application title
    titlePanel("Explore gene expression over the principal curve"),

    fluidPage(
      fluidRow(
        wellPanel(
          selectizeInput(
            inputId = 'GeneList', label = 'Select gene',
            choices = c('', sort(unique(unlist(lapply(lapply(InputList, "[[", "Expression"), colnames))))),
            selected = '',
            multiple = FALSE),
          lapply(1:nSpans, function(i){
            sliderInput(inputId = paste0("Span", i),
                        label = paste(InputList[[i]]$Name, "Span"), min = 0, max = 1, value = .25, step = .01)
            }
            )
        ),
        plotOutput("distPlot", height = PlotH)
      )
    )

  )

  # Define server logic
  server <- function(input, output, session) {

    # Plotting function ----------------------------------------

    output$distPlot <- renderPlot({

      gName <- input$GeneList

      print(paste("Plotting", gName))

      SpanVect <- rep(.25, length(InputList))

      for(i in 1:length(InputList)){
        SpanVect[i] <- input[[paste0("Span", i)]]
      }

      print(paste("Spans:", SpanVect))

      if(gName == ''){
        return(NULL)
      }

      plotList <- list()

      for(i in 1:length(InputList)){

        print(paste(i, 1))

        p1 <- PlotOnPseudotime(WorkStruct = InputList[[i]]$OrderedData,
                               Expression = InputList[[i]]$Expression,
                               Name = InputList[[i]]$Name,
                               gName = gName,
                               SpanVal = SpanVect[i],
                               CatOrder = NULL)

        if(gName %in% rownames(InputList[[i]]$Expression)){
          print(paste(i, "2a"))
          p2 <- PlotPG(X = InputList[[i]]$FinalExpMat,
                       TargetPG = InputList[[i]]$PGStruct,
                       GroupsLab = unlist(InputList[[i]]$Expression[rownames(InputList[[i]]$FinalExpMat),gName]),
                       p.alpha = .6,
                       Main = InputList[[i]]$Name,
                       PcToPlot = 1:2)
        } else {
          print(paste(i, "2b"))
          p2 <- PlotPG(X = InputList[[i]]$FinalExpMat,
                       TargetPG = InputList[[i]]$PGStruct,
                       GroupsLab = InputList[[i]]$FinalGroup,
                       p.alpha = .6,
                       Main = InputList[[i]]$Name,
                       PcToPlot = 1:2)
        }

        plotList[[length(plotList)+1]] <- p2[[1]]
        plotList[[length(plotList)+1]] <- p1

      }

      gridExtra::grid.arrange(grobs = plotList, ncol=2)

    })

  }


  # Run the application
  shinyApp(ui = ui, server = server)

}

