mfdmedian <- function(x,
                      type = "hdepth",
                      crossDepthsX = NULL,
                      depthOptions = NULL,
                      centerOption = "maxdepth") {

  ######
  # Check input.
  if (missing(x)) {
    stop("Input argument x is required.")
  }

  # Check x
  if (!is.array(x)) {
    stop("x must be a three dimensional array.")
  }
  if (length(dim(x)) != 3) {
    stop("x must be a three dimensional array.")
  }
  if (sum(is.nan(x)) != 0) {
    stop("x contains missing cases.")
  }
  t1 <- dim(x)[1]
  n1 <- dim(x)[2]
  p1 <- dim(x)[3]

  # Check type
  Indtype <- match(type, c("hdepth", "projdepth",
                           "sprojdepth", "dprojdepth", "sdepth"))[1]
  if (is.na(Indtype)) {
    stop("type should be one of hdepth, projdepth, sprojdepth, dprojdepth or sdepth")
  }
  if (Indtype == 5 & p1 > 2) {
    stop("Simplicial depth is only implemented for p<=2.")
  }

  # Check the crossDepths
  if (!is.null(crossDepthsX)) {
    if (!is.matrix(crossDepthsX)) {
      stop("crossDepthsX should be a matrix")
    }
    t2 <- ncol(crossDepthsX)
    n2 <- nrow(crossDepthsX)
    if (t1 != t2) {
      stop(paste("A different number of time points was specified ",
                 "between x and crossDepthsX."))
    }
    if (n1 != n2) {
      stop(paste("A different number of depths from the number of ",
                 "observations was specified."))
    }
    if (sum(colSums(crossDepthsX > 1)) != 0) {
      stop("The user supplied depths must take values in ]0,1].")
    }
    if (sum(colSums(crossDepthsX <= 0)) != 0) {
      stop("The user supplied depths must take values in ]0,1].")
    }
  }

  # check depthOptions
  if (is.null(depthOptions)) {
    depthOptions <- list()
  }
  if (!is.list(depthOptions)) {
    stop("depthOptions must be a list")
  }

  # check centerOption
  indCenter <- match(centerOption, c("maxdepth", "Huber"))[1]
  if (is.na(indCenter)) {
    stop("centerOption should be one of maxdepth or Huber.")
  }
  if (indCenter == 2 && Indtype != 2) {
    stop("Huber weights are only available for projection median.")
  }

  MFDmedian <- matrix(NaN, nrow = t1, ncol = p1)

  warningFlagFit <- 0
  warningIndFit <- c()

  options(warn = -1)
  for (j in 1:t1) {
    exactfit <- 0

    # R has standard dimension dropping, we need to be careful
    if (p1 == 1)  {
      xTimePoint <- matrix(x[j, , 1])
    }
    else {
      xTimePoint <- x[j, , , drop = TRUE]
    }

    # Find cross-sectional depth
    if (type == "hdepth") {
      temp <- try(do.call(hdepthmedian, c(list(x = xTimePoint), depthOptions)),
                  silent = TRUE)
      if (!is.list(temp)) {
        temp <- list()
      }
      if (!is.null(temp$median)) {
        MFDmedian[j, ] <- temp$median
      }
      else {
        exactfit <- 1
      }
    }
    else if (type == "projdepth") {
      if (is.null(crossDepthsX)) {
        temp <- projmedian(x = xTimePoint, options = depthOptions)
      } else {
        temp <- projmedian(x = xTimePoint,
                           projection.depths = crossDepthsX[, j])
      }
      if (!is.null(temp$max)) {
        if (indCenter == 1) {
          MFDmedian[j, ] <- temp$max
        }
        else {
          MFDmedian[j, ] <- temp$Huber
        }
      }
      else {
        exactfit <- 1
      }
    }
    else if (type == "sprojdepth") {
      if (is.null(crossDepthsX)) {
        temp <- sprojmedian(x = xTimePoint, options = depthOptions)
      } else {
        temp <- sprojmedian(x = xTimePoint,
                            sprojection.depths = crossDepthsX[, j])
      }
      if (!is.null(temp$max)) {
        if (indCenter == 1) {
          MFDmedian[j, ] <- temp$max
        }
      }
      else {
        exactfit <- 1
      }
    }
    else if (type == "dprojdepth") {
      if (is.null(crossDepthsX)) {
        temp <- dprojmedian(x = xTimePoint, options = depthOptions)
      } else {
        temp <- dprojmedian(x = xTimePoint,
                            dprojection.depths =  crossDepthsX[,j])
      }
      if (!is.null(temp$max)) {
        if (indCenter == 1) {
          MFDmedian[j, ] <- temp$max
        }
      }
      else {
        exactfit <- 1
      }
    }    
    else {
      temp <- sdepth(x = xTimePoint)
      if (!is.null(temp$depth)) {
        Ind <- which(temp$depth == max(temp$depth))
        MFDmedian[j, ] <- colMeans(matrix(xTimePoint[Ind, ], ncol = p1))
      }
      else {
        exactfit <- 1
      }
    }

    if (exactfit == 1) {
      warningFlagFit <- 1
      warningIndFit <- c(warningIndFit, j)
    }

  }
  options(warn = 0)

  Result <- list(MFDmedian = MFDmedian)
  class(Result) <- c("mrfDepth", "mfdmedian")
  

  # Handle all warnings
  if (warningFlagFit == 1) {
    warning(paste("Exact fits were detected for certain time points.",
                  "Their weights will be set to zero.",
                  "Check the returned results", call. = FALSE)
            )
    Result$IndFlagExactFit <- warningIndFit
  }

  return(Result)

}
