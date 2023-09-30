####    coef    ####

#' Coefficients method for (O)PLS models
#'
#' Coefficients of the (O)PLS(-DA) regression model
#'
#' @aliases coef.opls coef,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @return Numeric matrix of coefficients (number of rows equals the number of
#' variables, and the number of columns equals the number of responses)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' head(coef(sacurine.plsda))
#'
#' detach(sacurine)
#'
#' @rdname coef
#' @export
setMethod("coef", "opls",
          function(object) {
            return(object@coefficientMN)
          })


####    fitted    ####

#' Fitted method for 'opls' objects
#'
#' Returns predictions of the (O)PLS(-DA) model on the training dataset
#'
#' @aliases fitted.opls fitted,opls-method
#' @param object An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @return Predictions (either a vector, factor, or matrix depending on the y
#' response used for training the model)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' fitted(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname fitted
#' @export
setMethod("fitted", "opls",
          function(object) {
            
            if (!is.null(object@suppLs[["yPreMN"]])) {
              
              if (mode(object@suppLs[["yMCN"]]) == "character") {
                
                yPredMCN <- object@suppLs[[".char2numF"]](object@suppLs[["yPreMN"]],
                                                          c2nL = FALSE)
                
                if (is.vector(object@suppLs[["y"]])) {
                  fit <- c(yPredMCN)
                  names(fit) <- rownames(yPredMCN)
                } else if (is.factor(object@suppLs[["y"]])) {
                  fit <- c(yPredMCN)
                  names(fit) <- rownames(yPredMCN)
                  fit <- factor(fit, levels = levels(object@suppLs[["y"]]))
                } else if (is.matrix(object@suppLs[["y"]])) {
                  fit <- yPredMCN
                } else
                  stop() ## this case should not happen
                
              } else {
                
                yPredMCN <- object@suppLs[["yPreMN"]]
                
                if (is.vector(object@suppLs[["y"]])) {
                  fit <- c(yPredMCN)
                  names(fit) <- rownames(yPredMCN)
                } else if (is.matrix(object@suppLs[["y"]])) {
                  fit <- yPredMCN
                } else
                  stop() ## this case should not happen
                
              }
              
              return(fit)
              
            } else
              return(NULL)
            
          }) ## fitted



####    getLoadingMN    ####

#' @rdname getLoadingMN
#' @export
setMethod("getLoadingMN", "opls",
          function(object, orthoL = FALSE) {
            if(orthoL)
              return(object@orthoLoadingMN)
            else
              return(object@loadingMN)
          })




####    getPcaVarVn    ####

#' @rdname getPcaVarVn
#' @export
setMethod("getPcaVarVn", "opls",
          function(object) {
            return(object@pcaVarVn)
          })


####    getScoreMN    ####

#' @rdname getScoreMN
#' @export
setMethod("getScoreMN", "opls",
          function(object, orthoL = FALSE) {
            if (orthoL)
              return(object@orthoScoreMN)
            else
              return(object@scoreMN)
          })


####    getSubsetVi    ####

#' @rdname getSubsetVi
#' @export
setMethod("getSubsetVi", "opls",
          function(object) {
            return(object@subsetVi)
          })


####    getSummaryDF    ####
#'
#' @rdname getSummaryDF
#' @export
setMethod("getSummaryDF", "opls",
          function(object) {
            return(object@summaryDF)
          })


####    getVipVn    ####

#' @rdname getVipVn
#' @export
setMethod("getVipVn", "opls",
          function(object, orthoL = FALSE) {
            if(orthoL)
              return(object@orthoVipVn)
            else
              return(object@vipVn)
          })


####    getWeightMN    ####

#' @rdname getWeightMN
#' @export
setMethod("getWeightMN", "opls",
          function(object, orthoL = FALSE) {
            if(orthoL)
              return(object@orthoWeightMN)
            else
              return(object@weightMN)
          })


####    opls  (matrix)  ####

#' @rdname opls
#' @export
setMethod("opls", signature(x = "matrix"),
          function(x,
                   y = NULL,
                   predI = NA,
                   orthoI = 0,
                   
                   algoC = c("default", "nipals", "svd")[1],
                   crossvalI = 7,
                   log10L = FALSE,
                   permI = 20,
                   scaleC = c("none", "center", "pareto", "standard")[4],
                   subset = NULL,
                   
                   plotSubC = NA,                   
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink(info.txtC, append = TRUE)
            
            ## Checking arguments
            
            ## x -> xMN
            
            if (is.data.frame(x)) {
              if (!all(sapply(x, data.class) == "numeric")) {
                stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
              } else
                x <- as.matrix(x)
            } else if (is.matrix(x)) {
              if (mode(x) != "numeric")
                stop("'x' matrix must be of 'numeric' mode", call. = FALSE)
            } else
              stop("'x' must be either a data.frame or a matrix", call. = FALSE)
            
            if (any(apply(x, 2, function(colVn) all(is.na(colVn)))))
              stop("'x' contains columns with 'NA' only", call. = FALSE)
            
            if (any(is.infinite(x)))
              stop("'x' contains infinite values", call. = FALSE)
            
            if (any(is.nan(x)))
              stop("'x' contains NaN values", call. = FALSE)
            
            xMN <- x
            
            ## y -> yMCN
            
            yLevelVc <- NULL
            
            if (!is.null(y)) {
              
              if (is.vector(y)) {
                
                if (!(mode(y) %in% c("character", "numeric")))
                  stop("'y' vector must be of 'character' or 'numeric' type", call. = FALSE)
                
                if (length(y) != nrow(xMN))
                  stop("'y' vector length must be equal to the number of rows of 'x'", call. = FALSE)
                
                yMCN <- matrix(y, ncol = 1)
                
              } else if (is.factor(y)) {
                
                if (length(y) != nrow(xMN))
                  stop("'y' factor length must be equal to the number of rows of 'x'", call. = FALSE)
                
                yLevelVc <- levels(y)
                
                yMCN <- matrix(as.character(y), ncol = 1)
                
              } else if (is.matrix(y)) {
                
                if (!(mode(y) %in% c("character", "numeric")))
                  stop("'y' matrix must be of 'character' or 'numeric' type", call. = FALSE)
                
                if (nrow(y) != nrow(xMN))
                  stop("'x' and 'y' must have the same number of rows", call. = FALSE)
                
                if (ncol(y) > 1 && mode(y) != "numeric")
                  stop("Multiple response 'y' matrices must be of numeric type", call. = FALSE)
                
                yMCN <- y
                
              } else
                stop("'y' must be either a vector, a factor, or a matrix", call. = FALSE)
              
            } else
              yMCN <- NULL
            
            
            ## NA in Y only possible for multi-response regression (i.e., Y is a numeric matrix)
            
            if (!is.null(yMCN) &&
                ncol(yMCN) == 1 &&
                any(is.na(drop(yMCN))))
              stop("In case of single response modeling, 'y' must not contain missing values", call. = FALSE)
            
            if (!is.logical(log10L))
              stop("'log10L' must be a logical", call. = FALSE)
            
            if (permI < 0 || (permI - floor(permI)) > 1e-10)
              stop("'permI' must be an integer", call. = FALSE)
            
            if (permI > 0 && (is.null(yMCN) || ncol(yMCN) > 1)) {
              ## warning("Permutation testing available for single response (O)PLS(-DA) models only", call. = FALSE)
              permI <- 0
            }
            
            if (permI > 0 && !is.null(subset)) {
              permI <- 0
              warning("'permI' set to 0 because train/test partition is selected", call. = FALSE)
            }
            
            if (!(algoC %in% c('default', 'nipals', 'svd')))
              stop("'algoC' must be either 'default', 'nipals', or 'svd'", call. = FALSE)
            
            if (algoC == "default")
              algoC <- ifelse(is.null(yMCN) && !any(is.na(c(xMN))), "svd", "nipals")
            
            if (!is.null(yMCN) && algoC != "nipals")
              stop("'nipals' algorithm must be used for (O)PLS(-DA)", call. = FALSE)
            
            if ((is.na(orthoI) || orthoI > 0) && is.null(yMCN))
              stop("'y' response cannot be NULL for OPLS(-DA) modeling", call. = FALSE)
            
            if (!is.null(yMCN)) {
              if (is.na(orthoI) || orthoI > 0) {
                if (ncol(yMCN) > 1) {
                  stop("OPLS regression only available for a single 'y' response", call. = FALSE)
                } else if (mode(yMCN) == "character" && length(unique(drop(yMCN))) > 2)
                  stop("OPLS-DA only available for binary classification (use PLS-DA for multiple classes)", call. = FALSE)
              }
            }
            
            if (is.na(orthoI) || orthoI > 0)
              if (is.na(predI) || predI > 1) {
                predI <- 1
                warning("OPLS: number of predictive components ('predI' argument) set to 1", call. = FALSE)
              }
            
            if (!is.na(predI) && !is.na(orthoI) && ((predI + orthoI) > min(dim(xMN))))
              stop("The sum of 'predI' (", predI, ") and 'orthoI' (", orthoI, ") exceeds the minimum dimension of the 'x' data matrix (", min(dim(xMN)), ")" , call. = FALSE)
            
            if (!(length(scaleC) == 1 && scaleC %in% c('none', 'center', 'pareto', 'standard')))
              stop("'scaleC' must be either 'none', 'center', 'pareto', or 'standard'", call. = FALSE)
            
            if (!is.null(subset) && (is.null(yMCN) || ncol(yMCN) > 1))
              stop("train/test partition with 'subset' only available for (O)PLS(-DA) models of a single 'y' response", call. = FALSE)
            
            if (!is.null(subset) &&
                !(mode(subset) == 'character' && subset == 'odd') &&
                !all(subset %in% 1:nrow(xMN)))
              stop("'subset' must be either set to 'odd' or an integer vector of 'x' row numbers", call. = FALSE)
            
            if (crossvalI > nrow(xMN))
              stop("'crossvalI' must be less than the row number of 'x'", call. = FALSE)
            
            
            ## Constants
            
            epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16
            
            
            ## Character to numeric convertion function (for classification)
            
            if (!is.null(yMCN) && mode(yMCN) == "character") {
              
              if (!is.null(yLevelVc)) {
                claVc <- yLevelVc
              } else
                claVc <- sort(unique(drop(yMCN)))
              
              if (length(claVc) == 2) {
                ## binary response kept as a single vector for OPLS-DA computations
                .char2numF <- function(inpMCN,
                                       c2nL = TRUE) {
                  
                  if (c2nL) {
                    
                    outMCN <- inpMCN == claVc[2]
                    mode(outMCN) <- "numeric"
                    
                  } else {
                    
                    outMCN <- matrix(claVc[as.numeric(inpMCN > 0.5) + 1],
                                     ncol = 1,
                                     dimnames = dimnames(inpMCN))
                    
                  }
                  
                  return(outMCN)
                  
                }
                
              } else
                .char2numF <- function(inpMCN,
                                       c2nL = TRUE) {
                  
                  if (c2nL) {
                    
                    outMCN  <- t(sapply(drop(inpMCN),
                                        function(claC) as.numeric(claVc == claC)))
                    colnames(outMCN) <- claVc
                    
                    
                  } else {
                    
                    outMCN <- t(t(apply(inpMCN, 1,
                                        function(rowVn) claVc[which(rowVn == max(rowVn))[1]])))
                    colnames(outMCN) <- "y1"
                    
                  }
                  
                  return(outMCN)
                  
                }
              
            } else
              .char2numF <- NULL
            
            
            ####   Computations   ####
            
            
            ## rownames and colnames
            
            if (is.null(rownames(xMN)))
              if (!is.null(yMCN) && !is.null(rownames(yMCN))) {
                rownames(xMN) <- rownames(yMCN)
              } else
                rownames(xMN) <- paste0("s", 1:nrow(xMN))
            if (is.null(colnames(xMN)))
              colnames(xMN) <- paste0("x", 1:ncol(xMN))
            
            if (!is.null(yMCN)) {
              if (is.null(rownames(yMCN)))
                rownames(yMCN) <- rownames(xMN)
              if (is.null(colnames(yMCN)))
                colnames(yMCN) <- paste0("y", 1:ncol(yMCN))
            }
            
            ## Log10 transformation
            
            if (log10L)
              xMN <- .log10F(xMN)
            
            
            ## Test indices
            
            obsIniVi <- 1:nrow(xMN)
            
            if (!is.null(subset)) {
              subsetL <- TRUE
              if (length(subset) == 1 && subset == "odd") {
                
                if (mode(yMCN) == "numeric")
                  subsetVi <- seq(1, nrow(xMN), by = 2)
                else {
                  subsetVi <- integer()
                  for (claC in unique(drop(yMCN)))
                    subsetVi <- c(subsetVi,
                                  which(drop(yMCN) == claC)[seq(1, sum(drop(yMCN) == claC), by = 2)])
                  subsetVi <- sort(subsetVi)
                }
              } else
                subsetVi <- subset
              if (crossvalI > length(subsetVi))
                stop("'crossvalI' must be less than the number of samples in the subset", call. = FALSE)
            } else {
              subsetL <- FALSE
              subsetVi <- numeric()
            }
            
            
            ## Filtering out zero variance variables
            
            xVarIndLs <- list()
            xVarIndLs[[1]] <- 1:nrow(xMN)
            
            if (subsetL) {
              xVarIndLs[[1]] <- subsetVi
            } ## else if(!is.null(yMCN) && ncol(yMCN) == 1 && nrow(xMN) >= 2 * crossvalI)
            ##   for(cvkN in 1:crossvalI)
            ##       xVarIndLs <- c(xVarIndLs, list(setdiff(1:nrow(xMN), cvaOutLs[[cvkN]])))
            
            xVarVarLs <- lapply(xVarIndLs,
                                function(xVarVi) {
                                  apply(xMN[xVarVi, , drop = FALSE],
                                        2,
                                        function(colVn) var(colVn, na.rm = TRUE))
                                })
            
            xZeroVarVi <- integer()
            for (k in 1:length(xVarVarLs))
              xZeroVarVi <- union(xZeroVarVi, which(xVarVarLs[[k]] < epsN))
            
            if (length(xZeroVarVi) > 0) {
              names(xZeroVarVi) <- colnames(xMN)[xZeroVarVi]
              xMN <- xMN[, -xZeroVarVi, drop = FALSE]
              warning("The variance of the ",
                      length(xZeroVarVi),
                      " following variables is less than ",
                      signif(epsN, 2),
                      " in the full or partial (cross-validation) dataset: these variables will be removed:\n",
                      paste(names(xZeroVarVi), collapse = ", "),
                      call. = FALSE)
            }
            
            
            ## Core
            
            opl <- .coreF(xMN = xMN,
                          yMCN = yMCN,
                          orthoI = orthoI,
                          predI = predI,
                          scaleC = scaleC,
                          algoC = algoC,
                          crossvalI = crossvalI,
                          subsetL = subsetL,
                          subsetVi = subsetVi,
                          .char2numF = .char2numF,
                          info.txtC = info.txtC)
            
            if (is.null(y)) {          
              opl@suppLs["y"] <- list(NULL)
            } else
              opl@suppLs[["y"]] <- y
            
            if (is.null(opl@suppLs[["yMCN"]])) {
              opl@typeC <- "PCA"
            } else {
              if (ncol(opl@suppLs[["yMCN"]]) > 1 || mode(opl@suppLs[["yMCN"]]) == "numeric")
                opl@typeC <- "PLS"
              else
                opl@typeC <- "PLS-DA"
            }
            if (is.na(orthoI) || orthoI > 0)
              opl@typeC <- paste0("O", opl@typeC)
            
            opl@xZeroVarVi <- xZeroVarVi
            ## opl@suppLs[["yLevelVc"]] <- yLevelVc
            
            
            ## Permutation testing (Szymanska et al, 2012)
            
            if (permI > 0 && cumprod(dim(opl@summaryDF))[2] > 0) {
              
              modSumVc <- colnames(opl@summaryDF)
              
              permMN <- matrix(0,
                               nrow = 1 + permI,
                               ncol = length(modSumVc),
                               dimnames = list(NULL, modSumVc))
              
              perSimVn <- numeric(1 + permI)
              perSimVn[1] <- 1
              
              
              permMN[1, ] <- as.matrix(opl@summaryDF)
              
              for (k in 1:permI) {
                
                yVcn <- drop(opl@suppLs[["yMCN"]])
                if (!subsetL) {
                  yPerVcn <- sample(yVcn)
                } else {
                  yPerVcn <- numeric(nrow(xMN))
                  refVi <- opl@subsetVi
                  tesVi <- setdiff(1:nrow(xMN), refVi)
                  yPerVcn[refVi] <- sample(yVcn[refVi])
                  yPerVcn[tesVi] <- yVcn[tesVi]
                }
                yPerMCN <- matrix(yPerVcn, ncol = 1)
                
                perOpl <- .coreF(xMN = xMN,
                                 yMCN = yPerMCN,
                                 orthoI = opl@summaryDF[, "ort"],
                                 predI = opl@summaryDF[, "pre"],
                                 scaleC = scaleC,
                                 algoC = algoC,
                                 crossvalI = crossvalI,
                                 subsetL = subsetL,
                                 subsetVi = opl@subsetVi,
                                 .char2numF = .char2numF,
                                 info.txtC = info.txtC)
                
                permMN[1 + k, ] <- as.matrix(perOpl@summaryDF)
                
                perSimVn[1 + k] <- .similarityF(opl@suppLs[["yMCN"]], yPerMCN,
                                                .char2numF = .char2numF,
                                                charL = mode(opl@suppLs[["yMCN"]]) == "character")
                
              }
              
              permMN <- cbind(permMN, sim = perSimVn)
              
              perPvaVn <- c(pR2Y = (1 + length(which(permMN[-1, "R2Y(cum)"] >= permMN[1, "R2Y(cum)"]))) / (nrow(permMN) - 1),
                            pQ2 = (1 + length(which(permMN[-1, "Q2(cum)"] >= permMN[1, "Q2(cum)"]))) / (nrow(permMN) - 1))
              opl@summaryDF[, "pR2Y"] <- perPvaVn["pR2Y"]
              opl@summaryDF[, "pQ2"] <- perPvaVn["pQ2"]
              
              opl@suppLs[["permMN"]] <- permMN
              
            }
            
            ####   Numerical results   ####
            
            opl@descriptionMC <- rbind(samples = ifelse(!subsetL,
                                                        nrow(xMN),
                                                        length(subsetVi)),
                                       X_variables = ncol(xMN),
                                       near_zero_excluded_X_variables = length(opl@xZeroVarVi))
            
            totN <- length(c(xMN))
            nasN <- sum(is.na(c(xMN)))
            
            if (!is.null(opl@suppLs[["yMCN"]])) {
              
              opl@descriptionMC <- rbind(opl@descriptionMC,
                                         Y_variables = ncol(opl@suppLs[["yMCN"]]))
              totN <- totN + length(c(opl@suppLs[["yMCN"]]))
              nasN <- nasN + sum(is.na(c(opl@suppLs[["yMCN"]])))
              
            }
            
            opl@descriptionMC <- rbind(opl@descriptionMC,
                                       missing_values = paste0(nasN, " (", round(nasN / totN * 100), "%)"))
            
            
            ## Raw summary
            
            opl@suppLs[["topLoadI"]] <- 3
            
            if (ncol(xMN) > opl@suppLs[["topLoadI"]]) {
              xVarVn <- apply(xMN, 2, var)
              names(xVarVn) <- 1:length(xVarVn)
              xVarVn <- sort(xVarVn)
              xVarSorVin <- as.numeric(names(xVarVn[seq(1, length(xVarVn), length = opl@suppLs[["topLoadI"]])]))
              opl@suppLs[["xSubIncVarMN"]] <- xMN[, xVarSorVin, drop = FALSE]
            } else
              opl@suppLs[["xSubIncVarMN"]] <- xMN
            
            if (ncol(xMN) <= 100) {
              
              xCorMN <- cor(xMN, use = "pairwise.complete.obs")
              xCorMN[lower.tri(xCorMN, diag = TRUE)] <- 0
              
              if (ncol(xMN) > opl@suppLs[["topLoadI"]]) {
                
                xCorNexDF <- which(abs(xCorMN) >= sort(abs(xCorMN), decreasing = TRUE)[opl@suppLs[["topLoadI"]] + 1],
                                   arr.ind = TRUE)
                
                xCorDisMN <- matrix(0,
                                    nrow = nrow(xCorNexDF),
                                    ncol = nrow(xCorNexDF),
                                    dimnames = list(colnames(xMN)[xCorNexDF[, "row"]],
                                                    colnames(xMN)[xCorNexDF[, "col"]]))
                
                for (k in 1:nrow(xCorDisMN))
                  xCorDisMN[k, k] <- xCorMN[xCorNexDF[k, "row"], xCorNexDF[k, "col"]]
                
              } else
                xCorDisMN <- xCorMN
              
              opl@suppLs[["xCorMN"]] <- xCorDisMN
              
              rm(xCorDisMN)
              
            }
            
            ## Printing
            
            if (info.txtC != "none") {
              show(opl)
            }
            
            ## Plotting
            
            if (fig.pdfC != "none" && cumprod(dim(opl@summaryDF))[2] > 0)
              plot(opl, typeVc = "summary", plotSubC = plotSubC, fig.pdfC = fig.pdfC)
            
            ## Closing connection
            
            if (!(info.txtC %in% c("none", "interactive")))
              sink()
            
            ## Returning
            
            return(invisible(opl))
            
          })


#### opls (data.frame) ####

#' @rdname opls
#' @export
setMethod("opls", signature(x = "data.frame"),
          function(x,
                   y = NULL,
                   predI = NA,
                   orthoI = 0,
                   
                   algoC = c("default", "nipals", "svd")[1],
                   crossvalI = 7,
                   log10L = FALSE,
                   permI = 20,
                   scaleC = c("none", "center", "pareto", "standard")[4],
                   subset = NULL,
                   
                   plotSubC = NA,                   
                   fig.pdfC = c("none", "interactive", "myfile.pdf")[2],                   
                   info.txtC = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!all(sapply(x, data.class) == "numeric")) {
              stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
            } else
              x <- as.matrix(x)
            
            opl <- opls(x = x,
                        y = y,
                        predI = predI,
                        orthoI = orthoI,
                        
                        algoC = algoC,
                        crossvalI = crossvalI,
                        log10L = log10L,
                        permI = permI,
                        scaleC = scaleC,
                        subset = subset,
                        
                        plotSubC = plotSubC,                   
                        fig.pdfC = fig.pdfC,                   
                        info.txtC = info.txtC)
            
            return(invisible(opl))
            
          })


####    predict    ####

#' Predict method for (O)PLS models
#'
#' Returns predictions of the (O)PLS(-DA) model on a new dataset
#'
#' @aliases predict.opls predict,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param newdata Either a data frame or a matrix, containing numeric columns
#' only, with the same number of columns (variables) as the 'x' used for model
#' training with 'opls'.
#' @return Predictions (either a vector, factor, or matrix depending on the y
#' response used for training the model)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' predictorMN <- dataMatrix
#' responseFc <- sampleMetadata[, "gender"]
#'
#' sacurine.plsda <- opls(predictorMN,
#'                        responseFc,
#'                        subset = "odd")
#'
#' trainVi <- getSubsetVi(sacurine.plsda)
#'
#' table(responseFc[trainVi], fitted(sacurine.plsda))
#'
#' table(responseFc[-trainVi],
#'       predict(sacurine.plsda, predictorMN[-trainVi, ]))
#'
#' detach(sacurine)
#'
#' @rdname predict
#' @export
setMethod("predict", "opls",
          function(object, newdata) {
            
            if(object@typeC == "PCA")
              stop("Predictions currently available for (O)PLS(-DA) models only (not PCA)",
                   call. = FALSE)
            
            if(missing(newdata)) {
              
              return(fitted(object))
              
            } else {
              
              if(is.data.frame(newdata)) {
                if(!all(sapply(newdata, data.class) == "numeric")) {
                  stop("'newdata' data frame must contain numeric columns only", call. = FALSE)
                } else
                  newdata <- as.matrix(newdata)
              } else if(is.matrix(newdata)) {
                if(mode(newdata) != "numeric")
                  stop("'newdata' matrix must be of 'numeric' mode", call. = FALSE)
              } else
                stop("'newdata' must be either a data.frame or a matrix", call. = FALSE)
              
              if(ncol(newdata) != as.numeric(object@descriptionMC["X_variables", ])) {
                if(length(object@xZeroVarVi) == 0) {
                  stop("'newdata' number of variables is ",
                       ncol(newdata),
                       " whereas the number of variables used for model training was ",
                       as.numeric(object@descriptionMC["X_variables", ]),
                       ".",
                       call. = FALSE)
                } else if(ncol(newdata) - as.numeric(object@descriptionMC["X_variables", ]) ==
                          as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ])) {
                  warning(as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ]),
                          " near zero variance variables excluded during the model training will be removed from 'newdata'.",
                          call. = FALSE)
                  newdata <- newdata[, -object@xZeroVarVi, drop = FALSE]
                } else {
                  stop("'newdata' number of variables (",
                       ncol(newdata),
                       ") does not correspond to the number of initial variables (",
                       as.numeric(object@descriptionMC["X_variables", ]),
                       ") minus the number of near zero variance variables excluded during the training (",
                       as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ]),
                       ").",
                       call. = FALSE)
                }
              }
              
              xteMN <- scale(newdata, object@xMeanVn, object@xSdVn)
              
              if(object@summaryDF[, "ort"] > 0) {
                
                for(noN in 1:object@summaryDF[, "ort"]) {
                  if(object@suppLs[["naxL"]]) {
                    xtoMN <- matrix(0, nrow = nrow(xteMN), ncol = 1)
                    for(i in 1:nrow(xtoMN)) {
                      comVl <- complete.cases(xteMN[i, ])
                      xtoMN[i, ] <- crossprod(xteMN[i, comVl], object@orthoWeightMN[comVl, noN]) / drop(crossprod(object@orthoWeightMN[comVl, noN]))
                    }
                  } else
                    xtoMN <- xteMN %*% object@orthoWeightMN[, noN]
                  
                  xteMN <- xteMN - tcrossprod(xtoMN, object@orthoLoadingMN[, noN])
                }
                
              }
              
              if(object@suppLs[["naxL"]]) {
                yTesScaMN <- matrix(0, nrow = nrow(xteMN), ncol = ncol(object@coefficientMN),
                                    dimnames = list(rownames(xteMN), colnames(object@coefficientMN)))
                for(j in 1:ncol(yTesScaMN))
                  for(i in 1:nrow(yTesScaMN)) {
                    comVl <- complete.cases(xteMN[i, ])
                    yTesScaMN[i, j] <- crossprod(xteMN[i, comVl], object@coefficientMN[comVl, j])
                  }
              } else
                yTesScaMN <- xteMN %*% object@coefficientMN
              
              ## if(object@suppLs[["nayL"]])
              ##     yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]
              
              yTesMN <- scale(scale(yTesScaMN,
                                    FALSE,
                                    1 / object@ySdVn),
                              -object@yMeanVn,
                              FALSE)
              attr(yTesMN, "scaled:center") <- NULL
              attr(yTesMN, "scaled:scale") <- NULL
              
              if(is.factor(fitted(object))) {
                
                yTestMCN <- object@suppLs[[".char2numF"]](yTesMN,
                                                          c2nL = FALSE)
                predMCNFcVcn <- as.character(yTestMCN)
                names(predMCNFcVcn) <- rownames(newdata)
                predMCNFcVcn <- factor(predMCNFcVcn, levels = levels(object@suppLs[["y"]]))
                
              } else if(is.vector(fitted(object))) {
                
                if(is.character(fitted(object))) {
                  
                  yTestMCN <- object@suppLs[[".char2numF"]](yTesMN,
                                                            c2nL = FALSE)
                  predMCNFcVcn <- as.character(yTestMCN)
                  names(predMCNFcVcn) <- rownames(newdata)
                  
                } else {
                  
                  predMCNFcVcn <- as.numeric(yTesMN)
                  names(predMCNFcVcn) <- rownames(newdata)
                  
                }
                
              } else if (is.matrix(fitted(object))) {
                
                if (mode(fitted(object)) == "character") {
                  predMCNFcVcn  <- object@suppLs[[".char2numF"]](yTesMN,
                                                                 c2nL = FALSE)
                } else
                  predMCNFcVcn <- yTesMN
                
                rownames(predMCNFcVcn) <- rownames(newdata)
                
              }
              
              return(predMCNFcVcn)
              
            }
            
          })


####    print    ####

#' Print method for 'opls' objects
#'
#' Displays information about the dataset and the model.
#'
#' @aliases print.opls print,opls-method
#' @param x An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @return Invisible.
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' print(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname print
#' @export
setMethod("print", "opls",
          function(x) {
            
            cat("\n1) Data set:\n", sep = "")
            
            cat(x@descriptionMC["samples", ],
                " samples x ",
                x@descriptionMC["X_variables", ],
                " variables",
                ifelse(grepl("PLS", x@typeC),
                       paste0(" and ", ncol(x@suppLs[["yMCN"]]), " response", ifelse(ncol(x@suppLs[["yMCN"]]) > 1, "s", "")),
                       ""),
                "\n", sep = "")
            
            cat(x@descriptionMC["missing_values", ], " NAs\n", sep = "")
            
            cat(x@descriptionMC["near_zero_excluded_X_variables", ], " excluded variables (near zero variance)\n", sep = "")
            
            cat(x@suppLs[["scaleC"]], " x", ifelse(x@typeC != "PCA", " and y", ""), " scaling\n", sep = "")
            
            cat("Summary of the ", x@suppLs[["topLoadI"]], " increasing variance spaced raw variables:\n", sep = "")
            print(summary(x@suppLs[["xSubIncVarMN"]]))
            
            if (!is.null(x@suppLs[["xCorMN"]])) {
              cat("Correlations between the X-variables:\n")
              print(signif(x@suppLs[["xCorMN"]], 2))
              cat("\n", sep = "")
            }
            
            cat("\n2) Model: ", x@typeC, "\n", sep = "")
            
            cat("Correlations between variables and first 2 components:\n", sep = "")
            
            if (x@summaryDF[, "pre"] + x@summaryDF[, "ort"] < 2) {
              
              warning("A single component model has been selected by cross-validation", call. = FALSE)
              
              tCompMN <- x@scoreMN
              pCompMN <- x@loadingMN
              
            } else {
              
              if (x@summaryDF[, "ort"] > 0) {
                tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, 1])
                pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, 1])
                colnames(pCompMN) <- colnames(tCompMN) <- c("h1", "o1")
              } else {
                tCompMN <- x@scoreMN[, 1:2]
                pCompMN <- x@loadingMN[, 1:2]
              }
              
            }
            
            cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,
                             use = "pairwise.complete.obs")
            
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
            print(signif(topLoadMN, 2))
            
            message("")
            
            optDigN <- options()[["digits"]]
            options(digits = 3)
            print(x@modelDF)
            options(digits = optDigN)
            
          })


####    residuals    ####

#' Residuals method for (O)PLS models
#'
#' Returns the residuals from the (O)PLS(-DA) regression models
#'
#' @aliases residuals.opls residuals,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Numeric matrix or vector (same dimensions as the modeled y
#' response); if y is a character vector or a factor (in case of
#' classification), the residuals equal 0 (predicted class identical to the
#' true class) or 1 (prediction error)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.pls <- opls(dataMatrix,
#'                      sampleMetadata[, "age"])
#'
#' head(residuals(sacurine.pls))
#'
#' detach(sacurine)
#'
#' @rdname residuals
#' @export
setGeneric("residuals", getGeneric("residuals", package = "stats"))
setMethod("residuals", "opls",
          function(object, ...) {
            
            if (grepl("PLS", object@typeC)) {
              
              fit <- fitted(object)
              
              if (length(object@subsetVi) == 0) {
                trainVi <- 1:length(fit)
              } else
                trainVi <- object@subsetVi
              
              if (mode(object@suppLs[["yMCN"]]) == "numeric") {
                y <- object@suppLs[["y"]]
                if (is.matrix(y))
                  y <- y[trainVi, , drop = FALSE]
                else
                  y <- y[trainVi]
                return(y - fit)
              } else
                return(as.numeric(as.character(c(object@suppLs[["y"]])[trainVi]) != as.character(c(fit))))
              
            } else
              stop("'residuals' defined for (O)PLS(-DA) regression models only", call. = FALSE)
            
          }) ## residuals


####    show    ####

#' Show method for 'opls' objects
#'
#' Displays information about the dataset and the model.
#'
#' @aliases show.opls show,opls-method
#' @param object An S4 object of class \code{opls}, created by the \code{opls}
#' function.
#' @return Invisible.
#' @author Philippe Rinaudo and Etienne Thevenot (CEA)
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' show(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname show
#' @export
setMethod("show", "opls",
          function(object) {
            
            cat(object@typeC, "\n", sep = "")
            
            cat(object@descriptionMC["samples", ],
                " samples x ",
                object@descriptionMC["X_variables", ],
                " variables",
                ifelse(grepl("PLS", object@typeC),
                       paste0(" and ", ncol(object@suppLs[["yMCN"]]), " response", ifelse(ncol(object@suppLs[["yMCN"]]) > 1, "s", "")),
                       ""), "\n", sep = "")
            
            cat(object@suppLs[["scaleC"]], " scaling of predictors",
                ifelse(object@typeC == "PCA",
                       "",
                       paste0(" and ",
                              ifelse(mode(object@suppLs[["yMCN"]]) == "character" && object@suppLs[["scaleC"]] != "standard",
                                     "standard scaling of ",
                                     ""),
                              "response(s)")), "\n", sep = "")
            
            if (substr(object@descriptionMC["missing_values", ], 1, 1) != "0")
              cat(object@descriptionMC["missing_values", ], " NAs\n", sep = "")
            
            if (substr(object@descriptionMC["near_zero_excluded_X_variables", ], 1, 1) != "0")
              cat(object@descriptionMC["near_zero_excluded_X_variables", ],
                  " excluded variables (near zero variance)\n", sep = "")
            
            if (cumprod(dim(object@summaryDF))[2] > 0) {
              
              optDigN <- options()[["digits"]]
              options(digits = 3)
              print(object@summaryDF)
              options(digits = optDigN)
              
            } else {
              
              cat("Empty 'opls' object\n")
              
            }
            
          }) ## show


####    tested    ####
#'
#' @rdname tested
#' @export
setMethod("tested", "opls",
          function(object) {
            
            if(!is.null(object@suppLs[["yTesMN"]])) {
              
              if(mode(object@suppLs[["yMCN"]]) == "character") {
                
                yTestMCN <- object@suppLs[[".char2numF"]](object@suppLs[["yTesMN"]],
                                                          c2nL = FALSE)
                if(is.vector(object@suppLs[["y"]])) {
                  test <- c(yTestMCN)
                  names(test) <- rownames(yTestMCN)
                } else if(is.factor(object@suppLs[["y"]])) {
                  test <- c(yTestMCN)
                  names(test) <- rownames(yTestMCN)
                  test <- factor(test, levels = levels(object@suppLs[["y"]]))
                } else if(is.matrix(object@suppLs[["y"]])) {
                  test <- yTestMCN
                } else
                  stop() ## this case should not happen
                
              } else {
                
                yTestMCN <- object@suppLs[["yTesMN"]]
                
                if(is.vector(object@suppLs[["y"]])) {
                  test <- c(yTestMCN)
                  names(test) <- rownames(yTestMCN)
                } else if(is.matrix(object@suppLs[["y"]])) {
                  test <- yTestMCN
                } else
                  stop() ## this case should not happen
                
              }
              
              return(test)
              
            } else
              stop("Test results only available for (O)PLS(-DA) models", call. = FALSE)
            
            
          })