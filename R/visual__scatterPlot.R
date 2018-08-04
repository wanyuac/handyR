#' @title Making a scatter plot of x against y, colour coded by z.
#' @description Model: y ~ x.
#' @param y response variable of your interest to be plotted.
#' @param x explanatory variable of your interest to be plotted.
#' @param z the variable for generating a colour gradient. It is in the same order of both x and y.
#' @param x.bg and y.bg: background values of X and Y respectively to be contrasted with x and y. The function makes a scatter plot of x.bg and y.bg at first.
#' @param z.bg the variable for coding colours of (X.bg, Y.bg) with the same colour gradient defined by c(z, z.bg).
#' In practice, x, y and z are three columns of the same data frame.
#' log = c("both", "x", "y", "none"), determines whether to apply logarithmic transformation to x or y.
#' @param log.base the base for logarithmic transformation. The default transformtion takes log10.
#' @param replace.zero a small value used for replacing zero values before taking logarithms. Only takes effect if log != "none".
#' @param breaks.min and breaks.max determine the lower bound and upper bound for factorising the variable z.
#' @param breaks.n the number of breaks for dividing z into levels. The number of colour gradients is n-1.
#' @param col.grad a vector of characters defining boundaries of the colour gradient.
#' @param log.z whether to take the same logarithmic transformation for both z and z.bg.
#' @param log.z.base the base for the logarithmic transformation of z and z.bg.
#' @param output the name of the output PNG file.
#' @param w define the width of the output picture.
#' @param h define the height of the output picture.
#' @param title configure the title of the plot.
#' @param x.lab configure the label of X axis of the plot.
#' @param y.lab configure the label of Y axis of the plot.
#' @param x.extension and y.extension configure the extensions to the left/lower and right/upper sides of the plot.
#' @param margin a vector of integers to be passed to par(mar = margin).
#' For example, x.extension = c(-1, 1).
#' @param cex.main, cex.lab and cex.axis set the sizes of the title, labels and axes of the plot.
#' @param joint.col whether to colour all data points with the same colour gradient. Remember to check the break.min and breaks.max when you are using this option.
#' @param col.bg pch.bg, cex.bg configure the colour, style and size of the background data points (Y.bg ~ X.bg).
#' @param col.bg will be overridden if joint.col == TRUE.
#' @param pch and cex define the point style and size of data points of interests (Y ~ X).
#' @param show.range whethter to plot values from which an x is calculated. In this case, X is a summary statistic.
#' @param x.range a vector of characters separated by commas, containing every value that generating an X. For example, x.range = c("1,2,3", "2,2,5")
#' @param x.bg.range the same kind of vector as x.range for background X values.
#' @param range.pch and cex.range determine the type and size of points for values in x.range.
#' @param add.line determines whether to add a guideline to the plot.
#' @param line.a/b/h/v four arguments (a, b, h, v) to be passed to R's build-in function abline.
#' @param line.style the lty argument for abline(.)
#' @param line.col the colour of the guideline.
#' @param turnoff.dev not to call dev.off() if turnoff.dev = FALSE so that a user can add more features to the plot.
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Development history: 19 Nov 2015

scatterPlot <- function(x, y, z, x.bg = NA, y.bg = NA, z.bg = NA,
                        log = "both", log.base = 10, replace.zero = 1e-5,
                        breaks.min = 0, breaks.max = 1, breaks.n = 51,
                        col.grad = c("blue", "green", "red"),
                        log.z = FALSE, log.z.base = 10,
                        output = "plot.png", w = 800, h = 640,
                        title = "Y~X", x.lab = "X", y.lab = "Y",
                        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
                        x.extension = c(0, 0), y.extension = c(0, 0), margin = c(4.5, 4.5, 3, 2),
                        joint.col = FALSE, col.bg = "grey", pch.bg = 1, cex.bg = 1,
                        pch = 19, cex = 1,
                        show.range = FALSE, x.range = NA, x.bg.range = NA, range.pch = ".", cex.range = 1,
                        add.line = FALSE, line.a = NULL, line.b = NULL,
                        line.h = NULL, line.v = NULL, line.style = 3, line.col = "grey",
                        turnoff.dev = TRUE) {
  # If x.bg or y.bg is not provided, just replace x.bg and y.bg with x and y respectively.
  if (is.na(x.bg[1]) | is.na(y.bg[1])) {
    x.bg <- x
    y.bg <- y
    pch.bg <- pch
    cex.bg <- cex
  }

  # logarithmic transformation for X and Y
  if (log != "none") {
    if (log == "both") {
      x[x <= 0] <- replace.zero
      y[y <= 0] <- replace.zero
      x.bg[x.bg <= 0] <- replace.zero
      y.bg[y.bg <= 0] <- replace.zero
      x <- log(x, base = log.base)
      y <- log(y, base = log.base)
      x.bg <- log(x.bg, base = log.base)
      y.bg <- log(y.bg, base = log.base)
    } else if (log == "x") {
      x[x <= 0] <- replace.zero
      x.bg[x.bg <= 0] <- replace.zero
      x <- log(x, base = log.base)
      x.bg <- log(x.bg, base = log.base)
    } else {
      y[y <= 0] <- replace.zero
      y.bg[y.bg <= 0] <- replace.zero
      y <- log(y, base = log.base)
      y.bg <- log(y.bg, base = log.base)
    }
  }

  # logarithmic transformation for Z and Z.bg
  if (log.z) {
    z[z == 0] <- replace.zero  # replaces zero values before taking logarithms
    z.bg[z.bg == 0] <- replace.zero
    if (breaks.min == 0) {
      breaks.min <- replace.zero
    }
    if (breaks.max == 0) {
      breaks.max <- replace.zero
    }
    z <- log(z, base = log.z.base)
    z.bg <- log(z.bg, base = log.z.base)
    breaks.min == log(breaks.min, base = log.z.base)
    breaks.max == log(breaks.max, base = log.z.base)
  }

  # configure boundaries of the plot
  x.boundaries <- maxRange(x, x.bg)
  y.boundaries <- maxRange(y, y.bg)
  x.boundaries <- c(x.boundaries[1] + x.extension[1], x.boundaries[2] + x.extension[2])
  y.boundaries <- c(y.boundaries[1] + y.extension[1], y.boundaries[2] + y.extension[2])

  # set the colour gradients for (X, Y) and (X.bg, Y.bg)
  # use print(z.bins) to show the level of each one of z values
  # levels(z.bins): breaks.n breaking points and breaks.n - 1 bins
  # The function colorRampPalette(c("blue", "white", "red"), space = "rgb")(k) returns k colours ranging from blue to white to red.
  if (joint.col & !is.na(z.bg[1])) {
    z.joint <- c(z, z.bg)  # combine two vectors into a single one for making colour gradients
    z.bins <- cut(z.joint, breaks = seq(breaks.min, breaks.max, length.out = breaks.n), include.lowest = TRUE)
    palette <- colorRampPalette(col.grad, space = "rgb")(breaks.n - 1)
    z.n <- length(z)
    colours <- palette[z.bins[1 : z.n]]
    col.bg <- palette[z.bins[(z.n + 1) : length(z.bins)]]
  } else {
    z.bins <- cut(z, breaks = seq(breaks.min, breaks.max, length.out = breaks.n), include.lowest = TRUE)  # assign levels to values of z
    colours <- colorRampPalette(col.grad, space = "rgb")(breaks.n - 1)[z.bins]  # map colours to each z value based on its level
  }

  # create the plot and draw background data points
  png(filename = output, width = w, height = h)
  par(mar = margin)
  plot(x = x.bg, y = y.bg, xlim = x.boundaries, ylim = y.boundaries,
       main = title, xlab = x.lab, ylab = y.lab, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
       col = col.bg, pch = pch.bg, cex = cex.bg)

  # add ranges of X before adding (X, Y) to make the plot to show (X, Y) better
  if (show.range) {
    if (!is.na(x.bg.range[1])) {
      for (i in 1 : length(x.bg.range)) {
        x.bg.obs <- as.integer(unlist(strsplit(x.bg.range[i], ",")))
        y.bg.obs <- rep(y.bg[i], times = length(x.bg.obs))
        if (log != "none") {
          if (log == "both") {
            x.bg.obs[x.bg.obs <= 0] <- replace.zero
            y.bg.obs[y.bg.obs <= 0] <- replace.zero
            x.bg.obs <- log(x.bg.obs, base = log.base)
            y.bg.obs <- log(y.bg.obs, base = log.base)
          } else if (log == "x") {
            x.bg.obs[x.bg.obs <= 0] <- replace.zero
            x.bg.obs <- log(x.bg.obs, base = log.base)
          } else {
            y.bg.obs[y.bg.obs <= 0] <- replace.zero
            y.bg.obs <- log(y.bg.obs, base = log.base)
          }
        }
        points(x = x.bg.obs, y = y.bg.obs, col = col.bg, pch = range.pch, cex = cex.range)
      }
    }
    if (!is.na(x.range[1])) {
      for (i in 1 : length(x.range)) {
        x.obs <- as.integer(unlist(strsplit(x.range[i], ",")))
        y.obs <- rep(y[i], times = length(x.obs))
        if (log != "none") {
          if (log == "both") {
            x.obs[x.obs <= 0] <- replace.zero
            y.obs[y.obs <= 0] <- replace.zero
            x.obs <- log(x.obs, base = log.base)
            y.obs <- log(y.obs, base = log.base)
          } else if (log == "x") {
            x.obs[x.obs <= 0] <- replace.zero
            x.obs <- log(x.obs, base = log.base)
          } else {
            y.obs[y.obs <= 0] <- replace.zero
            y.obs <- log(y.obs, base = log.base)
          }
        }
        points(x = x.obs, y = y.obs, col = colours[i], pch = range.pch, cex = cex.range)
      }
    }
  }

  # add the guideline to the plot
  if (add.line) {
    abline(a = line.a, b = line.b, h = line.h, v = line.v, lty = line.style, col = line.col)
  }

  # Finally, add (X, Y) to the scatter plot.
  points(x = x, y = y, col = colours, pch = pch, cex = cex)

  if (turnoff.dev) {
    dev.off()
  }
}
