####    plot  (opls)  ####

#' Plot Method for (O)PLS(-DA)
#'
#' This function plots values based upon a model trained by \code{opls}.
#'
#' @aliases plot.opls plot,opls-method plot.oplsMultiDataSet plot,oplsMultiDataSet-method
#' @param x An S4 object of class \code{opls} or \code{oplsMultiDataSet},
#' created by the \code{opls} function.
#' @param y Currently not used
#' @param typeVc Character vector: the following plots are available:
#' 'correlation': Variable correlations with the components, 'outlier':
#' Observation diagnostics (score and orthogonal distances), 'overview': Model
#' overview showing R2Ycum and Q2cum (or 'Variance explained' for PCA),
#' 'permutation': Scatterplot of R2Y and Q2Y actual and simulated models after
#' random permutation of response values; 'predict-train' and 'predict-test':
#' Predicted vs Actual Y for reference and test sets (only if Y has a single
#' column), 'summary' [default]: 4-plot summary showing permutation, overview,
#' outlier, and x-score together, 'x-variance': Spread of raw variables
#' corresp. with min, median, and max variances, 'x-loading': X-loadings (the 6
#' of variables most contributing to loadings are colored in red to facilitate
#' interpretation), 'x-score': X-Scores, 'xy-score': XY-Scores, 'xy-weight':
#' XY-Weights
#' @param parAsColFcVn Optional factor character or numeric vector to be
#' converted into colors for the score plot; default is NA [ie colors will be
#' converted from 'y' in case of (O)PLS(-DA) or will be 'black' for PCA]
#' @param parCexN Numeric: amount by which plotting text should be magnified
#' relative to the default
#' @param parCompVi Integer vector of length 2: indices of the two components
#' to be displayed on the score plot (first two components by default)
#' @param parEllipsesL Should the Mahalanobis ellipses be drawn? If 'NA'
#' [default], ellipses are drawn when either a character parAsColVcn is
#' provided (PCA case), or when 'y' is a character factor ((O)PLS-DA cases).
#' @param parLabVc Optional character vector for the labels of observations on
#' the plot; default is NA [ie row names of 'x', if available, or indices of
#' 'x', otherwise, will be used]
#' @param parPaletteVc Optional character vector of colors to be used in the plots
#' @param parTitleL Should the titles of the plots be printed on the graphics
#' (default = TRUE); It may be convenient to set this argument to FALSE when
#' the user wishes to add specific titles a posteriori
#' @param parCexMetricN Numeric: magnification of the metrics at the bottom of
#' score plot (default -NA- is 1 in 1x1 and 0.7 in 2x2 display)
#' @param plotSubC Character: Graphic subtitle
#' @param fig.pdfC Character: File name with '.pdf' extension for the figure;
#' if 'interactive' (default), figures will be displayed interactively; if 'none',
#' no figure will be generated
#' @param info.txtC Character: File name with '.txt' extension for the printed
#' results (call to sink()'); if 'interactive' (default), messages will be
#' printed on the screen; if 'none', no verbose will be generated
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' for(typeC in c("correlation", "outlier", "overview",
#'                "permutation", "predict-train","predict-test",
#'                "summary", "x-loading", "x-score", "x-variance",
#'                "xy-score", "xy-weight")) {
#'
#'     print(typeC)
#'
#'     if(grepl("predict", typeC))
#'         subset <- "odd"
#'     else
#'         subset <- NULL
#'
#'     plsModel <- opls(dataMatrix, sampleMetadata[, "gender"],
#'                      predI = ifelse(typeC != "xy-weight", 1, 2),
#'                      orthoI = ifelse(typeC != "xy-weight", 1, 0),
#'                      permI = ifelse(typeC == "permutation", 10, 0),
#'                      subset = subset,
#'                      info.txtC = "none",
#'                      fig.pdfC = "none")
#'
#'     plot(plsModel, typeVc = typeC)
#'
#' }
#' 
#' sacPlsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#' plot(sacPlsda, parPaletteVc = c("green4", "magenta"))
#' 
#' detach(sacurine)
#'
#' @rdname plot
#' @export
setMethod("plot", signature(x = "opls"),
          function(x,
                   y,
                   typeVc = c("correlation",
                              "outlier",
                              "overview",
                              "permutation",
                              "predict-train",
                              "predict-test",
                              "summary",
                              "x-loading",
                              "x-score",
                              "x-variance",
                              "xy-score",
                              "xy-weight")[7],
                   parAsColFcVn = NA,
                   parCexN = 0.8,
                   parCompVi = c(1, 2),
                   parEllipsesL = NA,
                   parLabVc = NA,
                   parPaletteVc = NA,
                   parTitleL = TRUE,
                   parCexMetricN = NA,
                   plotSubC = "",
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (fig.pdfC == "none")
              stop("'fig.pdfC' cannot be set to 'none' in the 'plot' method.",
                   call. = FALSE)
            
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink(info.txtC, append = TRUE)
            
            if (length(typeVc) > 4) {
              warning("At most 4 graphics can be displayed simultaneously: here, the first 4 ones will be selected",
                      call. = FALSE)
              typeVc <- typeVc[1:4]
            }
            
            if ("summary" %in% typeVc) {
              if (!is.null(x@suppLs[["permMN"]]))
                typeVc <- c("overview",
                            "permutation",
                            "outlier",
                            "x-score")
              else
                typeVc <- c("overview",
                            "outlier",
                            "x-score",
                            "x-loading")
            }
            
            ## Checking arguments
             
            if (nrow(x@modelDF) < 1)
              stop("No model has been built and thus no plot can be displayed.", call. = FALSE)
            
            if (!all(typeVc %in% c('correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight')))
              stop("'typeVc' elements must be either 'correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight'", call. = FALSE)
            
            if ('predict-test' %in% typeVc && length(x@subsetVi) == 0)
              stop("For the 'predict-test' graphic to be generated, 'subset' must not be kept to NULL", call. = FALSE)
            
            if (!any(is.na(parLabVc))) {
              if (length(x@subsetVi) > 0 && length(parLabVc) != nrow(x@suppLs[["yMCN"]])) {
                stop("When 'subset' is not NULL, 'parLabVc' vector length must be equal to the number of train + test samples (here: ", nrow(x@suppLs[["yMCN"]]), ").", call. = FALSE)
              } else if (length(parLabVc) != nrow(x@scoreMN))
                stop("'parLabVc' vector length must be equal to the number of 'x' rows", call. = FALSE)
              if (mode(parLabVc) != "character")
                stop("'parLabVc' must be of 'character' type", call. = FALSE)
            }
            
            # eset <- getEset(x)
            # 
            # if (is.na(plotSubC) && cumprod(dim(Biobase::exprs(eset)))[2] > 0)
            #   plotSubC <- Biobase::experimentData(eset)@title
            # if (nchar(plotSubC) > 32)
            #   plotSubC <- paste0(substr(plotSubC, 1, 32), ".")
            #   
            # if (!is.na(plotPhenoDataC)) {
            #   if (cumprod(dim(Biobase::exprs(eset)))[2] < 1)
            #     stop("'plotPhenoDataC' can be used only for models computed on ExpressionSet objects")
            #   if (!is.character(plotPhenoDataC))
            #     stop("'plotPhenoDataC' must be a character when the 'plot' method is applied to an 'opls' instance")
            #   pdataDF <- Biobase::pData(eset)
            #   if (!(plotPhenoDataC %in% colnames(pdataDF))) {
            #     stop("'plotPhenoDataC' must be the name of a column of the sampleMetadata slot of the 'ExpressionSet' instance")
            #   } else
            #     parAsColFcVn <- pdataDF[, plotPhenoDataC]
            # }
            
            # if (!any(is.na(parAsColFcVn))) {
            if (!all(is.na(parAsColFcVn))) {
              if (length(x@subsetVi) > 0 && length(parAsColFcVn) != nrow(x@suppLs[["yMCN"]])) {
                stop("When 'subset' is not NULL, 'parAsColFcVn' vector length must be equal to the number of train + test samples (here: ", nrow(x@suppLs[["yMCN"]]), ").", call. = FALSE)
              } else if (length(parAsColFcVn) != nrow(x@scoreMN))
                stop("'parAsColFcVn' vector length must be equal to the number of 'x' rows", call. = FALSE)
                if (!(mode(parAsColFcVn) %in% c("character", "numeric")))
                stop("'parAsColFcVn' must be of 'character' or 'numeric' type", call. = FALSE)
              if (is.character(parAsColFcVn)) {
                parAsColFcVn <- factor(parAsColFcVn)
                # warning("Character 'parAsColFcVn' set to a factor", call. = FALSE)
              }
            }
            
            if (is.null(x@suppLs[["permMN"]]) && 'permutation' %in% typeVc)
              stop("'permI' must be > 0 for 'permutation' graphic to be plotted", call. = FALSE)
            
            if (x@summaryDF[, "ort"] > 0)
              if (parCompVi[1] != 1) {
                parCompVi[1] <- 1
                warning("OPLS: first component to display ('parCompVi' first value) set to 1", call. = FALSE)
              }
            
            if ("xy-weight" %in% typeVc &&
                substr(x@typeC, 1, 3) != "PLS")
              ## (is.null(yMCN) || is.na(x@summaryDF[, "ort"]) || x@summaryDF[, "ort"] > 0))
              stop("'xy-weight graphic can be displayed only for PLS(-DA) models", call. = FALSE)
            
            if (any(grepl('predict', typeVc)))
              if (is.null(x@suppLs[["yMCN"]]) ||
                  ncol(x@suppLs[["yMCN"]]) > 1 ||
                  (mode(x@suppLs[["yMCN"]]) == "character" && length(unique(drop(x@suppLs[["yMCN"]]))) > 2))
                ## if(any(grepl('predict', typeVc)) && is.matrix(x@fitted"]]) && ncol(x@fitted"]]) > 1)
                ## if(any(grepl('predict', typeVc)) && (is.null(yMCN) || ncol(yMCN) != 1))
                stop("'predict' graphics available for single response regression or binary classification only", call. = FALSE)
            
            if (is.na(parEllipsesL)) {
              if ((all(is.na(parAsColFcVn)) && grepl("-DA$", x@typeC)) ||
                  (!all(is.na(parAsColFcVn)) && is.factor(parAsColFcVn))) {
                parEllipsesL <- TRUE
              } else
                parEllipsesL <- FALSE
              ## if((x@typeC == "PCA" && !all(is.na(parAsColFcVn)) && is.factor(parAsColFcVn)) || ## PCA case
              ##    grepl("-DA$", x@typeC)) { ## (O)PLS-DA cases
              ##     parEllipsesL <- TRUE
              ## } else
              ##     parEllipsesL <- FALSE
            } else if (parEllipsesL && !grepl("-DA$", x@typeC) && (all(is.na(parAsColFcVn)) || !is.factor(parAsColFcVn)))
              stop("Ellipses can be plotted for PCA (or PLS regression) only if the 'parAsColFcVn' is a factor",
                   call. = FALSE)
            
            if (x@summaryDF[, "pre"] + x@summaryDF[, "ort"] < 2) {
              
              if (!all(typeVc %in% c("permutation", "overview"))) {
                warning("Single component model: only 'overview' and 'permutation' (in case of single response (O)PLS(-DA)) plots available", call. = FALSE)
                typeVc <- "overview"
                if (!is.null(x@suppLs[["permMN"]]))
                  typeVc <- c(typeVc, "permutation")
              }
              
              tCompMN <- x@scoreMN
              pCompMN <- x@loadingMN
              
            } else {
              
              if (x@summaryDF[, "ort"] > 0) {
                if (parCompVi[2] > x@summaryDF[, "ort"] + 1)
                  stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
                tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, parCompVi[2] - 1])
                pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, parCompVi[2] - 1])
                colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
              } else {
                if (max(parCompVi) > x@summaryDF[, "pre"])
                  stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
                tCompMN <- x@scoreMN[, parCompVi, drop = FALSE]
                pCompMN <- x@loadingMN[, parCompVi, drop = FALSE]
              }
              
            }
            
            ## if(ncol(tCompMN) > 1) {
            
            ##     mahInvCovMN <- solve(cov(tCompMN))
            
            ##     pcaResMN <- cbind(sdsVn = apply(tCompMN,
            ##                           1,
            ##                           function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
            ##                       odsVn = apply(x@suppLs[["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
            ##                           1,
            ##                           function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))
            
            ## } else
            ##     pcaResMN <- NULL
            
            cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,
                             use = "pairwise.complete.obs")
            
            if (!is.null(x@suppLs[["yModelMN"]]))
              cytCompMN <- cor(x@suppLs[["yModelMN"]], tCompMN, use = "pairwise.complete.obs")
            
            
            if (x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]])) {
              
              pexVi <- integer(x@suppLs[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values
              
              for (k in 1:ncol(pCompMN)) {
                
                pkVn <-  pCompMN[, k]
                
                pexVi[1:(2 * x@suppLs[["topLoadI"]]) + 2 * x@suppLs[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x@suppLs[["topLoadI"]]],
                                                                                                  rev(order(pkVn, decreasing = TRUE)[1:x@suppLs[["topLoadI"]]]))
                
              }
              
            } else
              pexVi <- 1:ncol(x@suppLs[["xModelMN"]])
            
            pxtCompMN <- cbind(pCompMN,
                               cxtCompMN)
            
            if (ncol(pCompMN) == 1) {
              colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
            } else
              colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])
            
            topLoadMN <- pxtCompMN
            
            topLoadMN <- topLoadMN[pexVi, , drop = FALSE]
            
            if (x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]]) &&
                ncol(pCompMN) > 1) {
              
              topLoadMN[(2 * x@suppLs[["topLoadI"]] + 1):(4 * x@suppLs[["topLoadI"]]), c(1, 3)] <- NA
              topLoadMN[1:(2 * x@suppLs[["topLoadI"]]), c(2, 4)] <- NA
              
            }
            
            if (!any(is.na(parPaletteVc))) {
              
              if (mode(parPaletteVc) != "character")
                stop("'parPaletteVc' should be a vector of mode 'character'",
                     call. = FALSE)
              
              palChkVl <- sapply(parPaletteVc,
                                 function(colC) {
                                   tryCatch(is.matrix(grDevices::col2rgb(colC)), 
                                            error = function(e) FALSE)
                                   
                                 }) ## as proposed by Josh O'Brien on stackoverflow
              
              if (any(is.na(palChkVl)) || any(!palChkVl))
                stop("The following elements from 'parPaletteVc' could not be interpreted as colors:\n",
                     paste(names(palChkVl)[is.na(palChkVl) | !palChkVl], collapse = ", "),
                     call. = FALSE)
              
              if (any(is.na(palChkVl)))
                stop("The following element(s) from 'parPaletteVc' could not be interpreted as colors:\n",
                     paste(names(palChkVl)[is.na(palChkVl)], collapse = ", "),
                     call. = FALSE)
              
            }
            
            ## Observation and variable names and colors
            
            ## obsLabVc
            
            if (!any(is.na(parLabVc))) {
              obsLabVc <- parLabVc
            } else if (!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
              obsLabVc <- rownames(x@suppLs[["yMCN"]])
            } else {## PCA
              if (!is.null(rownames(tCompMN))) {
                obsLabVc <- rownames(tCompMN)
              } else
                obsLabVc <- as.character(1:nrow(tCompMN))
            }
            
            if (length(x@subsetVi) > 0) {
              ## (O)PLS(-DA) models of a single 'y' response
              tesLabVc <- obsLabVc[-x@subsetVi]
              obsLabVc <- obsLabVc[x@subsetVi]
            } else
              tesLabVc <- ""
            
            ## obsColVc

            if (!all(is.na(parAsColFcVn))) {
              obsColVc <- .plotColorF(as.vector(parAsColFcVn), parPaletteVc)[["colVc"]]
              obsLegVc <- as.vector(parAsColFcVn)
            } else if (!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
              obsColVc <- .plotColorF(c(x@suppLs[["yMCN"]]), parPaletteVc)[["colVc"]]
              obsLegVc <- c(x@suppLs[["yMCN"]])
            } else {## PCA
              obsColVc <- rep("black", nrow(tCompMN))
              obsLegVc <- NULL
            }            
            # if (!any(is.na(parAsColFcVn))) {
            #   obsColVc <- .plotColorF(as.vector(parAsColFcVn), parPaletteVc)[["colVc"]]
            #   obsLegVc <- as.vector(parAsColFcVn)
            # } else if (!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
            #   obsColVc <- .plotColorF(c(x@suppLs[["yMCN"]]), parPaletteVc)[["colVc"]]
            #   obsLegVc <- c(x@suppLs[["yMCN"]])
            # } else {## PCA
            #   obsColVc <- rep("black", nrow(tCompMN))
            #   obsLegVc <- NULL
            # }
            
            if (length(x@subsetVi) > 0) {
              ## (O)PLS(-DA) models of a single 'y' response
              tesColVc <- obsColVc[-x@subsetVi]
              obsColVc <- obsColVc[x@subsetVi]
              if (!is.null(obsLegVc)) {
                tesLegVc <- obsLegVc[-x@subsetVi]
                obsLegVc <- obsLegVc[x@subsetVi]
              }
            }
            
            ## Plotting parameters
            
            if (fig.pdfC != "interactive")
              pdf(fig.pdfC)
            
            if (length(typeVc) > 1) {
              opar <- par(mfrow = c(2, 2),
                          mar = c(4.6, 4.1, 2.6, 1.6),
                          font = 2,
                          font.axis = 2,
                          font.lab = 2,
                          lwd = 2,
                          pch = 18)
              layL <- TRUE
            } else {
              opar <- par(mar = c(5.1, 4.1, 4.1, 2.1),
                          font = 2,
                          font.axis = 2,
                          font.lab = 2,
                          lwd = 2,
                          pch = 18)
              layL <- FALSE
            }
            
            ## Graph
            
            for (ploC in typeVc) {
              if (length(typeVc) == 1 ||
                  ploC == "overview") {
                parTitleSetC <- plotSubC
              } else
                parTitleSetC <- ""
              
              .plotF(ploC,
                     opl = x,
                     obsColVc = obsColVc,
                     obsLabVc = obsLabVc,
                     obsLegVc = obsLegVc,
                     layL = layL,
                     parCexN = parCexN,
                     parCexMetN = parCexMetricN,
                     parEllipsesL = parEllipsesL,
                     parPaletteVc = parPaletteVc,
                     parTitleL = parTitleL,
                     parCompVi = parCompVi,
                     typeVc = typeVc,
                     tCompMN = tCompMN,
                     pCompMN = pCompMN,
                     cxtCompMN = cxtCompMN,
                     cytCompMN = cytCompMN,
                     topLoadMN = topLoadMN,
                     pexVi = pexVi,
                     tesColVc = tesColVc,
                     tesLabVc = tesLabVc,
                     tesLegVc = tesLegVc,
                     parTitleSetC = parTitleSetC)
            }
            
            par(opar)
            
            if (fig.pdfC != "interactive")
              dev.off()
            
            ## Closing connection
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
          })

.plotF <- function(ploC,
                   opl,
                   obsColVc,
                   obsLabVc,
                   obsLegVc,
                   layL,
                   parCexN,
                   parCexMetN,
                   parEllipsesL,
                   parPaletteVc,
                   parTitleL,
                   parCompVi,
                   typeVc,
                   tCompMN,
                   pCompMN,
                   cxtCompMN,
                   cytCompMN,
                   topLoadMN,
                   pexVi,
                   tesColVc,
                   tesLabVc,
                   tesLegVc,
                   parTitleSetC) {
  
  ploPclF <- function() {
    
    xLimVn <- NULL
    yLimVn <- NULL
    
    ploColVc <- "black"
    
    if (ploC == "correlation") {
      
      maiC <- "Variable correlations"
      if (parTitleSetC != "")
        maiC <- paste0(maiC, "\n", parTitleSetC)
      
      xLabC <- paste("with t",
                     parCompVi[1],
                     sep = "")
      
      yLabC <- paste("with t",
                     parCompVi[2],
                     sep = "")
      
      if (opl@summaryDF[, "ort"] > 0)
        yLabC <- paste("with tOrtho",
                       parCompVi[2] - 1,
                       sep = "")
      
      yLimVn <- xLimVn <- c(-1, 1)
      
      ploMN <- cxtCompMN
      
      if (opl@typeC != "PCA")
        ploMN <- rbind(ploMN,
                       cytCompMN)
      
    } else if (substr(ploC, 1, 7) == "predict") {
      
      maiC <- paste0("Predicted vs Actual",
                     paste0(" (",
                            unlist(strsplit(ploC, "-"))[2],
                            ")"))
      if (parTitleSetC != "")
        maiC <- paste0(maiC, "\n", parTitleSetC)
      
      xLabC <- "predicted"
      yLabC <- "actual"
      
      if (grepl("train", ploC)) {
        ploNamVc <- obsLabVc
        ploColVc <- obsColVc
      } else {
        ploNamVc <- tesLabVc
        ploColVc <- tesColVc
      }
      
      ypMN <- eval(parse(text = paste("opl@suppLs[['y", switch(unlist(strsplit(ploC, "-"))[2], train = "Pre", test = "Tes"), "MN']]", sep = ""))) ## predicted
      
      if (length(opl@subsetVi) == 0) {
        yaMCN <- opl@suppLs[["yMCN"]] ## actual
      } else {
        if (grepl("train", ploC))
          yaMCN <- opl@suppLs[["yMCN"]][opl@subsetVi, , drop = FALSE]
        else
          yaMCN <- opl@suppLs[["yMCN"]][-opl@subsetVi, , drop = FALSE]
      }
      
      if (mode(opl@suppLs[["yMCN"]]) == "character") { ## binary only
        ypMN <- ypMN[, 1, drop = FALSE]
        yaMN <- opl@suppLs[[".char2numF"]](yaMCN)[, 1, drop = FALSE]
      } else
        yaMN <- yaMCN
      
      ploMN <- cbind(ypMN,
                     yaMN) ## to be modified (when ncol(yPreMCN) > 1)
      
    } else if (ploC == "x-loading") {
      
      maiC <- "Loadings"
      if (parTitleSetC != "")
        maiC <- paste0(maiC, "\n", parTitleSetC)
      
      xLabC <- paste0("p",
                      parCompVi[1],
                      " (",
                      round(opl@modelDF[parCompVi[1], "R2X"] * 100),
                      "%)")
      
      yLabC <- paste0("p",
                      parCompVi[2],
                      " (",
                      round(opl@modelDF[parCompVi[2], "R2X"] * 100),
                      "%)")
      
      ploMN <- pCompMN
      
      if (!is.null(opl@suppLs[["yMCN"]]) && opl@summaryDF[, "ort"] > 0)
        yLabC <- paste0("pOrtho",
                        parCompVi[2] - 1,
                        " (",
                        round(opl@modelDF[parCompVi[2] - 1, "R2X"] * 100),
                        "%)")
      
      
    } else if (ploC == "x-score") {
      
      maiC <- paste0("Scores (", opl@typeC, ")")
      if (parTitleSetC != "")
        maiC <- paste0(maiC, "\n", parTitleSetC)
      
      xLabC <- paste0("t",
                      parCompVi[1],
                      " (",
                      round(opl@modelDF[parCompVi[1], "R2X"] * 100),
                      "%)")
      
      yLabC <- paste0("t",
                      parCompVi[2],
                      " (",
                      round(opl@modelDF[parCompVi[2], "R2X"] * 100),
                      "%)")
      
      ploMN <- tCompMN
      
      if (grepl("^OPLS", opl@typeC))
        yLabC <- paste0("to", parCompVi[2] - 1)
      
      xLimVn <- c(-1, 1) * max(sqrt(var(ploMN[, 1]) * hotFisN), max(abs(ploMN[, 1])))
      yLimVn <- c(-1, 1) *max(sqrt(var(ploMN[, 2]) * hotFisN), max(abs(ploMN[, 2])))
      
      ploColVc <- obsColVc
      
    } else if (ploC == "xy-score") {
      
      maiC <- "XY-Scores"
      if (parTitleSetC != "")
        maiC <- paste0(maiC, "\n", parTitleSetC)
      
      xLabC <- paste("t", parCompVi[1], sep = "")
      yLabC <- paste("u/c", parCompVi[1], sep = "")
      
      ploMN <- cbind(opl@scoreMN[, parCompVi[1]], opl@uMN[, parCompVi[1]] / opl@cMN[parCompVi[1]])
      
      ploColVc <- obsColVc
      
    } else if (ploC == "xy-weight") {
      
      maiC <- "Weights"
      if (parTitleSetC != "")
        maiC <- paste0(maiC, "\n", parTitleSetC)
      xLabC <- paste0("w*c", parCompVi[1])
      yLabC <- paste0("w*c", parCompVi[2])
      
      ploMN <- rbind(opl@weightStarMN[, parCompVi],
                     opl@cMN[, parCompVi])
      
      pchVn <- rep(17, times = nrow(ploMN))
      ploColVc <- rep("grey", times = nrow(ploMN))
      
      pchVn[(nrow(opl@weightStarMN) + 1):nrow(ploMN)] <- 15
      ploColVc[(nrow(opl@weightStarMN) + 1):nrow(ploMN)] <- "black"
      
    }
    
    
    if (is.null(xLimVn))
      xLimVn <- range(ploMN[, 1])
    if (is.null(yLimVn))
      yLimVn <- range(ploMN[, 2])
    
    plot(ploMN,
         main = ifelse(parTitleL, maiC, ""),
         type = "n",
         xlab = xLabC,
         ylab = yLabC,
         xlim = xLimVn,
         ylim = yLimVn)
    
    abline(v = axTicks(1),
           col = "grey")
    
    abline(h = axTicks(2),
           col = "grey")
    
    abline(v = 0)
    abline(h = 0)
    
    if (ploC == "correlation") {
      
      lines(cos(radVn),
            sin(radVn))
      
      corPexVi <- pexVi
      corPchVn <- rep(18, nrow(cxtCompMN))
      corNamVc <- rownames(cxtCompMN)
      ## corPchVn <- rep(18, ncol(opl@suppLs[["xModelMN"]]))
      ## corNamVc <- colnames(opl@suppLs[["xModelMN"]])
      
      
      if(opl@typeC != "PCA") {
        corPexVi <- c(corPexVi, (nrow(cxtCompMN) + 1):nrow(ploMN))
        corPchVn <- c(corPchVn, rep(15, nrow(cytCompMN)))
        corNamVc <- c(corNamVc, rownames(cytCompMN))
      }
      
      points(ploMN,
             pch = corPchVn)
      
      points(ploMN[corPexVi, ],
             pch = corPchVn[corPexVi],
             col = "red")
      
      text(ploMN[corPexVi, ],
           cex = parCexN,
           col = "red",
           labels = corNamVc[corPexVi],
           pos = 3)
      
      if (opl@typeC != "PCA" && length(typeVc) == 1)
        legend("topleft",
               pch = c(18, 15),
               legend = c("X vars", "Y vars"))
      
    } else if (substr(ploC, 1, 7) == "predict") {
      
      abline(0, 1)
      
      text(ploMN[, 1:2],
           cex = parCexN,
           col = ploColVc,
           labels = ploNamVc)
      
      if (!is.null(obsLegVc))
        if (ploC == "predict-train") {
          .plotLegendF(obsLegVc,
                       ploMN,
                       paletteVc = parPaletteVc)
          
        } else
          .plotLegendF(tesLegVc,
                       ploMN,
                       paletteVc = parPaletteVc)
      
    } else if (ploC == "x-loading") {
      
      points(ploMN,
             col = "grey",
             pch = 18)
      
      points(ploMN[pexVi, ],
             pch = 18,
             col = "black")
      
      ## pexLabVc <- colnames(opl@suppLs[["xModelMN"]])[pexVi]
      pexLabVc <- rownames(opl@loadingMN)[pexVi]
      pexLabVc[duplicated(pexLabVc)] <- ""
      
      text(ploMN[pexVi, ],
           cex = parCexN,
           col = "black",
           labels = pexLabVc,
           pos = rep(c(4, 2, 3, 1), each = opl@suppLs[["topLoadI"]]))
      
    } else if (ploC == "x-score") {
      
      lines(sqrt(var(ploMN[, 1]) * hotFisN) * cos(radVn),
            sqrt(var(ploMN[, 2]) * hotFisN) * sin(radVn))
      ## Tenenhaus98, p87
      
      if (!is.null(obsLegVc))
        .plotLegendF(obsLegVc,
                     ploMN,
                     paletteVc = parPaletteVc)
      
      text(ploMN,
           cex = parCexN,
           col = ploColVc,
           labels = obsLabVc)
      
      pu1N <- par("usr")[1]
      pu2N <- par("usr")[2]
      
      if (is.na(parCexMetN))
        parCexMetN <- ifelse(layL, 0.7, 1)
      
      mtext(paste("R2X", round(opl@summaryDF[, "R2X(cum)"], 3), sep = "\n"),
            at = pu1N * ifelse(layL, 1.35, 1.1),
            cex = parCexMetN,
            font = 1,
            line = 3,
            side = 1)
      
      
      if (parEllipsesL) {
        par(lwd = 2)
        for (colC in unique(ploColVc)) {
          ploColVl <- ploColVc == colC
          if (sum(ploColVl) > 1)
            .plotEllipseF(ploMN[ploColVl, , drop = FALSE],
                          colC = colC)
        }
      }
      
      if (!is.null(opl@suppLs[["yMCN"]])) {
        
        mtext(paste("R2Y", round(opl@summaryDF[, "R2Y(cum)"], 3), sep = "\n"),
              at = pu1N * ifelse(layL, 1, 0.8),
              cex = parCexMetN,
              font = 1,
              line = 3,
              side = 1)
        
        mtext(paste("Q2Y", round(opl@summaryDF[, "Q2(cum)"], 3), sep = "\n"),
              at = pu1N * ifelse(layL, 0.6, 0.4),
              cex = parCexMetN,
              font = 1,
              line = 3,
              side = 1)
        
        mtext(paste("RMSEE", round(opl@summaryDF[, "RMSEE"], 3), sep = "\n"),
              at =  -pu1N * ifelse(layL, 0.6, 0.4),
              cex = parCexMetN,
              font = 1,
              line = 3,
              side = 1)
        
        mtext(paste("pre", opl@summaryDF[, "pre"], sep = "\n"),
              at = -pu1N * ifelse(layL, 0.92, 0.7),
              cex = parCexMetN,
              font = 1,
              line = 3,
              side = 1)
        
        if (opl@summaryDF[, "ort"] > 0)
          mtext(paste("ort", opl@summaryDF[, "ort"], sep = "\n"),
                at = -pu1N * ifelse(layL, 1.1, 0.9),
                cex = parCexMetN,
                font = 1,
                line = 3,
                side = 1)
        
      }
      
    } else if (ploC == "xy-score") {
      
      abline(0, 1)
      
      if (!is.null(obsLegVc))
        .plotLegendF(obsLegVc,
                     ploMN,
                     paletteVc = parPaletteVc)
      
      text(ploMN,
           cex = parCexN,
           col = ploColVc,
           labels = obsLabVc)
      
    } else if (ploC == "xy-weight") {
      
      text(ploMN[, 1:2],
           cex = parCexN,
           col = ploColVc,
           labels = c(rownames(opl@weightStarMN), rownames(opl@cMN)))
      
      if (!layL)
        legend("topleft",
               col = c("grey", "black"),
               legend = c("X", "Y"),
               text.col = c("grey", "black"))
      
    }
    
  } ## end of ploPclF()
  
  
  if (is.null(tCompMN) && ploC %in% c("correlation",
                                      "outlier",
                                      "x-loading",
                                      "x-score",
                                      "xy-weight"))
    warning("No ", ploC, " plotting", call. = FALSE)
  
  ## Hotteling's T2 (Tenenhaus98, p86)
  
  if (!is.null(tCompMN))
    hotFisN <- (nrow(tCompMN) - 1) * 2 * (nrow(tCompMN)^2 - 1) / (nrow(tCompMN) * nrow(tCompMN) * (nrow(tCompMN) - 2)) * qf(0.95, 2, nrow(tCompMN) - 2)
  
  
  radVn <- seq(0, 2 * pi, length.out = 100)
  
  
  if (ploC == "outlier") {
    
    ## Observation diagnostics
    ## see Hubert2005 p66
    
    mahInvCovMN <- solve(cov(tCompMN))
    
    pcaResMN <- cbind(sdsVn = apply(tCompMN,
                                    1,
                                    function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
                      odsVn = apply(opl@suppLs[["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
                                    1,
                                    function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))
    
    pcaResThrVn <- c(sqrt(qchisq(0.975, 2)),
                     (mean(pcaResMN[, 2]^(2/3)) + sd(pcaResMN[, 2]^(2/3)) * qnorm(0.975))^(3/2))
    
    pcaResExtVi <- union(which(pcaResMN[, 1] > pcaResThrVn[1]),
                         which(pcaResMN[, 2] > pcaResThrVn[2]))
    
    maiC <- "Observation diagnostics"
    if (parTitleSetC != "")
      maiC <- paste0(maiC, "\n", parTitleSetC)
    
    plot(pcaResMN,
         main = maiC,
         type = "n",
         xlab = "Score distance (SD)",
         xlim = c(0, max(pcaResMN[, 1]) * 1.1),
         xpd = TRUE,
         ylab = "Orthogonal distance (OD)",
         ylim = c(0, max(pcaResMN[, 2]) * 1.1))
    abline(v = pcaResThrVn[1],
           lty = "dashed")
    abline(h = pcaResThrVn[2],
           lty = "dashed")
    
    if (length(pcaResExtVi)) {
      
      points(pcaResMN[-pcaResExtVi, , drop = FALSE],
             col = obsColVc[-pcaResExtVi],
             pch = 18)
      text(pcaResMN[pcaResExtVi, , drop = FALSE],
           cex = parCexN,
           col = obsColVc[pcaResExtVi],
           labels = obsLabVc[pcaResExtVi])
      
    } else
      points(pcaResMN,
             col = obsColVc,
             pch = 18)
    
  } ## outlier
  
  
  if (ploC == "overview") {
    
    if (opl@typeC == "PCA") {
      
      barplot(opl@modelDF[, "R2X"] * 100,
              main = ifelse(parTitleSetC != "",
                            parTitleSetC,
                            "Explained variance"),
              names.arg = rownames(opl@modelDF),
              xlab = "PC",
              ylab = "% of total variance")
      
      
    } else {
      
      if (opl@summaryDF[, "ort"] == 0) {
        modBarDF <- opl@modelDF
      } else
        modBarDF <- opl@modelDF[!(rownames(opl@modelDF) %in% c("rot", "sum")), ]
      
      barplot(rbind(modBarDF[, "R2Y(cum)"],
                    modBarDF[, "Q2(cum)"]),
              beside = TRUE,
              main = ifelse(parTitleSetC != "",
                            parTitleSetC,
                            "Model overview"),
              names.arg = rownames(modBarDF),
              xlab = "")
      
      axis(2, lwd=2, lwd.ticks=1)
      
      abline(h = 0.5,
             col = "gray")
      
      barplot(rbind(modBarDF[, "R2Y(cum)"],
                    modBarDF[, "Q2(cum)"]),
              add = TRUE,
              beside = TRUE,
              col = c("grey", "black"))
      
      text(1.5,
           0,
           adj = c(0, 0.5),
           col = "white",
           srt = 90,
           labels = " R2Y") ## R2Ycum
      
      text(2.5,
           0,
           adj = c(0, 0.5),
           col = "white",
           labels = " Q2Y", ## Q2cum
           srt = 90)
      
    }
    
  } ## overview
  
  
  if (ploC == "permutation") {
    
    # perDiaVc <- c("R2Y(cum)", "Q2(cum)")
    
    maiC <- paste0("pR2Y = ",
                   opl@summaryDF[, "pR2Y"],
                   ", pQ2 = ",
                   opl@summaryDF[, "pQ2"])
    if (parTitleSetC != "")
      maiC <- paste0(maiC, "\n", parTitleSetC)
    
    plot(c(min(opl@suppLs[["permMN"]][, "sim"]), 1),
         c(min(opl@suppLs[["permMN"]][, c("R2Y(cum)", "Q2(cum)")]), 1),
         main = maiC,
         type = "n",
         xlab = expression(Similarity(bold(y), bold(y[perm]))),
         ylab = "")
    
    points(opl@suppLs[["permMN"]][, "sim"], opl@suppLs[["permMN"]][, "Q2(cum)"],
           col = "black",
           pch = 18)
    abline(h = opl@suppLs[["permMN"]][1, "Q2(cum)"],
           col = "black")
    points(opl@suppLs[["permMN"]][, "sim"], opl@suppLs[["permMN"]][, "R2Y(cum)"],
           col = "grey",
           pch = 18)
    abline(h = opl@suppLs[["permMN"]][1, "R2Y(cum)"],
           col = "grey")
    .plotLegendF(c("R2Y", "Q2Y"),
                 "bottomright",
                 colVc = c("grey", "black"))
    
  } ## permutation
  
  
  if (ploC == "x-variance") {
    
    par(las = 2)
    
    maiC <- "X variances (min, med, max)"
    if (parTitleSetC != "")
      maiC <- paste0(maiC, "\n", parTitleSetC)
    
    boxplot(opl@suppLs[["xSubIncVarMN"]],
            main = maiC,
            names = rep("", 3),
            xaxt = "n",
            yaxt = "n")
    
    axis(at = axTicks(2),
         side = 2)
    
    mtext(substr(colnames(opl@suppLs[["xSubIncVarMN"]]), 1, 9),
          cex = ifelse(layL, 0.7, 0.8),
          at = 1:3,
          line = 0.2,
          side = 1)
    
    par(las = 0)
    
  } ## "x-variance"
  
  
  if (ploC %in% c("correlation",
                  "predict-test",
                  "predict-train",
                  "x-loading",
                  "x-score",
                  "xy-score",
                  "xy-weight"))
    ploPclF()
  
} ## .plotF

## Transforms a character or numeric vector into colors
.plotColorF <- function(namVcn, palVc = NA) {
  
  if (any(is.na(palVc))) {
    ## 16 color palette without 'gray'
    palVc <- c("blue", "red", "green3", "cyan", "magenta", "#FF7F00", "#6A3D9A", "#B15928", "aquamarine4", "yellow4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#FFFF99")
  }
  
  if (is.null(namVcn) || all(is.na(namVcn))) {
    
    if (!is.null(namVcn)) {
      
      dev.new()
      
      palNamVc <- paste0(1:length(palVc),
                         "_",
                         palVc)
      
      pie(rep(1, length(palVc)),
          col = palVc,
          labels = palNamVc)
      
      print(matrix(palNamVc, ncol = 1))
      
    }
    
    return(palVc)
    
  } else {
    
    if (is.character(namVcn)) {
      
      namFcn <- factor(namVcn)
      
      if (length(levels(namFcn)) <= length(palVc)) {
        scaVc <- palVc[1:length(levels(namFcn))]
      } else
        scaVc <- c(palVc,
                   rep("gray",
                       length(levels(namFcn)) - length(palVc)))
      
      names(scaVc) <- levels(namFcn)
      
      colVc <- scaVc[unlist(sapply(namVcn,
                                   function(scaleC) {
                                     if (is.na(scaleC))
                                       return(NA)
                                     else
                                       which(levels(namFcn) == scaleC)
                                   }))]
      
    } else if (is.numeric(namVcn)) {
      
      scaVc <- rev(rainbow(100, end = 4/6))
      if (length(namVcn) > 1) {
        colVc <- scaVc[round((namVcn - min(namVcn, na.rm = TRUE)) / diff(range(namVcn, na.rm = TRUE)) * 99) + 1]
      } else
        colVc <- rep("black", length(namVcn))
      
    } else
      stop("'namVcn' argument must be a vector of either character or numeric mode", call. = FALSE)
    
    colVc[is.na(colVc)] <- "black"
    names(colVc) <- namVcn
    
  }
  
  return(list(colVc = colVc,
              scaVc = scaVc))
  
}


## Draws Mahalanobis ellipse
.plotEllipseF <- function(xMN,
                          colC = NULL,
                          sxyMN = NULL) {
  ## Adapted from the 'drawMahal' function of the 'chemometrics' package
  ## by P. Filzmoser and K. Varmuza
  
  if (ncol(xMN) != 2)
    stop("Matrix must have two columns", call. = FALSE)
  
  xMN <- xMN[!apply(xMN, 1, function(x) any(is.na(x))), , drop = FALSE]
  
  radVn <- seq(0, 2 * pi, length.out = 100)
  
  csqN <- qchisq(0.95, 2) ## ncol(xMN) == 2
  
  xMeaVn <- colMeans(xMN)
  ##        t1        t2
  ## 1.1771851 0.5661031
  
  xCovMN <- sxyMN
  
  if (is.null(xCovMN))
    xCovMN <- cov(xMN)
  ##            t1         t2
  ## t1  1.8079514 -0.9768156
  ## t2 -0.9768156  1.0201432
  
  xCovSvdLs <- svd(xCovMN, nv = 0)
  ## $ d: num [1:2] 2.467 0.361
  ## $ u: num [1:2, 1:2] -0.829 0.559 0.559 0.829
  ## $ v: NULL
  
  if (!is.null(colC)) {
    
    mahMN <- matrix(1, nrow = length(radVn), ncol = 1) %*% xMeaVn + cbind(cos(radVn), sin(radVn)) %*% diag(sqrt(xCovSvdLs[["d"]] * csqN)) %*% t(xCovSvdLs[["u"]])
    lines(mahMN,
          col = colC)
    
  } else {
    
    zerVarVin <- which(xCovSvdLs[["d"]] < .Machine$double.eps)
    
    if (length(zerVarVin))
      stop("Covariance matrix cannot be inverted because of ", length(zerVarVin), " zero eigen values\n", call. = FALSE)
    else
      sxyInvMN <- xCovSvdLs[["u"]] %*% diag(1 / xCovSvdLs[["d"]]) %*% t(xCovSvdLs[["u"]])
    
    invisible(sxyInvMN)
    
  }
  
} ## end of .plotEllipseF()


## Plots the figure legend
.plotLegendF <- function(namOrLegVcn,
                         locCMN = "topright",
                         txtCexN = 0.7,
                         colVc = NULL,
                         paletteVc = NA) {
  ## Note:
  ##  locCMN: either a character indicating the corner of the plot where the legend is to be plotted or the numeric matrix of point coordinates for the legLocF function below to find the corner where there is most space
  
  ## Determining the location (corner) for the legend
  
  legLocF <- function(thrN = 0.2) {
    
    lefN <- par("usr")[1] + thrN * diff(par("usr")[1:2])
    rigN <- par("usr")[2] - thrN * diff(par("usr")[1:2])
    topN <- par("usr")[4] - thrN * diff(par("usr")[3:4])
    botN <- par("usr")[3] + thrN * diff(par("usr")[3:4])
    
    locVl <- c(all(ploMN[, 1] > lefN |
                     ploMN[, 2] < topN),
               all(ploMN[, 1] < rigN |
                     ploMN[, 2] < topN),
               all(ploMN[, 1] > lefN |
                     ploMN[, 2] > botN),
               all(ploMN[, 1] < rigN |
                     ploMN[, 2] > botN))
    names(locVl) <- c("topleft", "topright", "bottomleft", "bottomright")
    
    return(locVl)
    
  }
  
  stopifnot(is.character(locCMN) || is.matrix(locCMN))
  
  if(is.matrix(locCMN)) {
    
    ploMN <- locCMN
    
    thrVn <- seq(0, 0.25, by = 0.05)
    locSumVn <- sapply(thrVn, function(thrN) sum(legLocF(thrN = thrN)))
    
    if(sum(locSumVn) > 0)
      locC <- names(which(legLocF(thrVn[max(which(locSumVn > 0))]))[1])
    else
      locC <- "topleft"
    
  } else
    locC <- locCMN
  
  
  ## Determining the color scale
  
  if(!is.null(colVc)) {
    scaVc <- colVc
    names(scaVc) <- namOrLegVcn
  } else
    scaVc <- .plotColorF(namOrLegVcn, paletteVc)[["scaVc"]]
  
  legTypC <- ifelse(is.character(namOrLegVcn), "cha", "num")
  
  ## Plotting the legend
  
  dpx <- diff(par("usr")[1:2])
  dpy <- diff(par("usr")[3:4])
  pu1 <- par("usr")[1]
  pu2 <- par("usr")[2]
  pu3 <- par("usr")[3]
  pu4 <- par("usr")[4]
  
  if(locC == "topright") {
    
    xLefN <- pu2 - 0.05 * dpx
    xRigN <- pu2 - 0.02 * dpx
    yBotN <- pu4 - 0.22 * dpy
    yTopN <- pu4 - 0.02 * dpy
    
  } else if(locC == "topleft") {
    
    xLefN <- pu1 + 0.02 * dpx
    xRigN <- pu1 + 0.05 * dpx
    yBotN <- pu4 - 0.22 * dpy
    yTopN <- pu4 - 0.02 * dpy
    
  } else if(locC == "bottomleft") {
    
    xLefN <- pu1 + 0.02 * dpx
    xRigN <- pu1 + 0.05 * dpx
    yBotN <- pu3 + 0.02 * dpy
    yTopN <- pu3 + 0.22 * dpy
    
  } else if(locC == "bottomright") {
    
    xLefN <- pu2 - 0.05 * dpx
    xRigN <- pu2 - 0.02 * dpx
    yBotN <- pu3 + 0.02 * dpy
    yTopN <- pu3 + 0.22 * dpy
    
  }
  
  yVn <- seq(yBotN,
             yTopN,
             length = length(scaVc) + 1)
  
  yBotVn <- yVn[1:(length(yVn) - 1)]
  yTopVn <- yVn[2:length(yVn)]
  
  rect(xleft = xLefN,
       ybottom = yBotVn,
       xright = xRigN,
       ytop = yTopVn,
       col = scaVc,
       border = NA)
  
  xLegN <- ifelse(grepl("left", locC), xRigN, xLefN)
  xAdjN <- ifelse(grepl("left", locC), 0, 1)
  
  if(legTypC == "cha") {
    text(xLegN,
         seq(yBotN + (yTopN - yBotN) / (2 * length(scaVc)),
             yTopN - (yTopN - yBotN) / (2 * length(scaVc)),
             length = length(scaVc)),
         adj = c(xAdjN, 0.5),
         cex = txtCexN,
         col = scaVc,
         labels = names(scaVc))
  } else
    text(xLegN,
         seq(yBotN,
             yTopN,
             length = 5),
         adj = c(xAdjN, 0.5),
         cex = txtCexN,
         labels = signif(seq(min(namOrLegVcn, na.rm = TRUE), max(namOrLegVcn, na.rm = TRUE), length = 5), 2))
  
}
