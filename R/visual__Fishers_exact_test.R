#' @title Simulating and plotting probabilities in a Fisher's exact test
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# Editions: 4-19 Oct 2015

#========== DEFINE THE FUNCTION ==========
plotFishersP <- function(r1, c1, n, alpha, panels = 3, prefix = "FishersProbability_", font.size = 1.5, w = 1200, h = 1000) {
  # n: sample size
  # r1: n(A), r2: n(B), and they should be smaller than or equal to n.
  # As a constrain by the contingency table, r1 + c1 <= n in the simulation.
  cat("n(A): ", r1, "\n", sep = "")
  cat("n(B): ", c1, "\n", sep = "")
  cat("sample size: ", n, "\n", sep = "")

  test.n <- min(r1, c1)
  accum.p <- exact.p <- OR <- numeric(test.n + 1)  # numeric(test.n + 1) includes a == 0
  denominator <- choose(n, c1)
  for (a in 0 : test.n) {  # You cannot get "a" greater than min(r1, c1)
    exact.p[a + 1] <- choose(r1, a) * choose(n - r1, c1 - a) / denominator
    m <- matrix(data = c(a, c1 - a, r1 - a, n + a - r1 - c1), nrow = 2, ncol = 2)
    rownames(m) <- colnames(m) <- c("present", "absent")
    t <- fisher.test(m, alternative = "two.sided")
    accum.p[a + 1] <- t$p.value
    OR[a + 1] <- t$estimate  # the conditional MLE OR
    cat("a = ", a, ", OR = ", OR[a + 1], ", P = ", accum.p[a + 1], "\n", sep = "")
  }
  x.series <- seq(from = 0, to = test.n, by = 1)

  png(filename = paste(prefix, paste(r1, c1, n, sep = "_"), ".png", sep = ""), width = w, height = h)
  if (panels == 3) {
    par(mfrow = c(2, 2), mar = c(4, 4.5, 4, 4))  # a two-by-two grid for plotting

    # PLOT EXACT P VALUES
    plot(x = x.series, y = exact.p, xlab = "a", ylab = "p", main = "Individual probabilities",
         cex.lab = font.size, cex.axis = font.size, cex.main = font.size + 0.5)
    #abline(v = r1, col = "red", lty = 3)  # incidence of A
    #abline(v = c1, col = "green", lty = 3)  # incidence of B
    #axis(side = 1, at = r1)  # add a tick at the X axis
    #axis(side = 1, at = c1)
    max.index <- which.max(exact.p) - 1  # find out the index of the maximum value
    abline(v = max.index, col = "green", lty = 3)
    axis(side = 1, at = max.index, cex.axis = font.size)
  } else {
    par(mfrow = c(1, 2), mar = c(4, 4.5, 4, 4))  # only a single row and two columns for plotting
  }

  # PLOT ORs
  plot(x = x.series, y = OR, xlab = "a", ylab = "Odds ratio", main = "Odds ratios",
       cex.lab = font.size, cex.axis = font.size, cex.main = font.size + 0.5)
  OR.one <- which.min(abs(OR - 1)) - 1  # find the one cloesest to OR = 1
  abline(h = 1, col = "red", lty = 3)  # to distinguish positive and negative associations
  axis(side = 4, at = 1, cex.axis = font.size)
  abline(v = OR.one, col = "green", lty = 3)
  axis(side = 1, at = OR.one, cex.axis = font.size)
  cat("Number of negative associations: ", sum(OR < 1), "\n", sep = "")

  # PLOT ACCUMULATED TWO-WAY P VALUES
  plot(x = x.series, y = accum.p, xlab = "a", ylab = "Two-sided P value", main = "Two-sided P values",
       cex.lab = font.size, cex.axis = font.size, cex.main = font.size + 0.5)
  abline(h = alpha, col = "red", lty = 3)  # the threshold for P values
  axis(side = 4, at = alpha, cex.axis = font.size)
  break.point <- which.min(abs(accum.p - alpha)) - 1  # the point closest to the thresold of P values
  abline(v = break.point, col = "green", lty = 3)
  axis(side = 1, at = break.point, cex.axis = font.size)

  dev.off()

  return(list(exact.p, accum.p))
}

#========== EXECUTION ==========
#n <- 300
#p <- plotFishersP(r1 = 80, c1 = 200, n = n, alpha = 0.05, panels = 3, prefix = "FishersProbability_3panels_", font.size = 2)  # incidences of A, B, and the sample size
#p <- plotFishersP(r1 = 80, c1 = 200, n = n, alpha = 0.05, panels = 2, prefix = "FishersProbability_2panels_", font.size = 2)
#p <- plotFishersP(r1 = 3, c1 = 5, n = n, alpha = 0.05, panels = 2, prefix = "FishersProbability_2panels_", font.size = 2)
#p <- plotFishersP(r1 = 5, c1 = 200, n = n, alpha = 0.05, panels = 2, prefix = "FishersProbability_2panels_", font.size = 2)
#p <- plotFishersP(r1 = 2, c1 = 3, n = n, alpha = 0.05, panels = 2, prefix = "FishersProbability_2panels_", font.size = 2)  # singleton
#p <- plotFishersP(r1 = 1, c1 = 1, n = n, alpha = 0.05, panels = 2, prefix = "FishersProbability_2panels_", font.size = 2)  # singleton
