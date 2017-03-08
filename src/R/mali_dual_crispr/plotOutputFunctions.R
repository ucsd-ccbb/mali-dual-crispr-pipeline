plot_fit <- function(timepointsList, x1, ac1, fc1, ab1, x2, ac2, fc2, ab2, nt, pA_pB, minfc = 0.10) {
  nx <- nrow(x1)
  maxt <- timepointsList[nt] + 3
  plot_lambda1 <- rep(0, maxt)
  plot_lambda2 <- rep(0, maxt)
  
  # TODO: Do these 1, 2 represent hardcoding of # of expected replicates?  If so, need to refactor function
  good1 <- t(t(x1) > ab1)
  good2 <- t(t(x2) > ab2)
  allbad <-
    apply(good1, 1, sum) < 2 &
    apply(good2, 1, sum) < 2 #in this case I have nothing to use in either experiment
  
  for (i in 1:maxt)
    plot_lambda1[i] <- -log2(sum(2 ^ (ac1 + fc1 * i)))
  for (i in 1:maxt)
    plot_lambda2[i] <- -log2(sum(2 ^ (ac2 + fc2 * i)))
  
  rge <- range(c(x1, x2))
  for (i in 1:nx) {
    if (allbad[i])
      next
    if (abs(fc1[i]) < minfc & abs(fc2[i]) < minfc)
      next
    pch <- rep(1, nt)
    pch[good1[i, ]] <- 16 #circle
    plot(
      timepointsList,
      x1[i, ],
      ylim = rge,
      pch = pch,
      cex = 1.2,
      main = pA_pB[i],
      xlab = "time",
      ylab = expression(paste(log[2], " relative frequency"))
    )
    yfit <- ac1[i] + fc1[i] * (1:maxt) + plot_lambda1
    lines(1:maxt, yfit)
    
    pch <- rep(1, nt)
    pch[good2[i, ]] <- 16
    points(timepointsList,
           x2[i, ],
           pch = pch,
           cex = 1.2,
           col = "blue")
    yfit <- ac2[i] + fc2[i] * (1:maxt) + plot_lambda2
    lines(1:maxt, yfit, col = "blue")
  }
}

plotOverlappingHist <-
  function(a,
           b,
           colors = c("gray70", "gray20", "gray50"),
           breaks = NULL,
           xlim = NULL,
           ylim = NULL,
           xlab = NULL,
           ylab = NULL,
           main = NULL) {
    ahist = NULL
    bhist = NULL
    if (!(is.null(breaks))) {
      ahist = hist(a, breaks = breaks, plot = FALSE)
      bhist = hist(b, breaks = breaks, plot = FALSE)
    } else {
      ahist = hist(a, plot = FALSE)
      bhist = hist(b, plot = FALSE)
      dist = ahist$breaks[2] - ahist$breaks[1]
      breaks = seq(min(ahist$breaks, bhist$breaks),
                   max(ahist$breaks, bhist$breaks),
                   dist)
      ahist = hist(a, breaks = breaks, plot = FALSE)
      bhist = hist(b, breaks = breaks, plot = FALSE)
    }
    
    if (is.null(xlim)) {
      xlim = c(min(ahist$breaks, bhist$breaks),
               max(ahist$breaks, bhist$breaks))
    }
    
    if (is.null(ylim)) {
      ylim = c(0, max(ahist$counts, bhist$counts))
    }
    
    overlap = ahist
    for (i in 1:length(overlap$counts)) {
      if (ahist$counts[i] > 0 & bhist$counts[i] > 0) {
        overlap$counts[i] = min(ahist$counts[i], bhist$counts[i])
      } else {
        overlap$counts[i] = 0
      }
    }
    
    plot(
      ahist,
      xlim = xlim,
      ylim = ylim,
      col = colors[1],
      border = colors[1],
      xlab = xlab,
      ylab = ylab,
      main = main
    )
    plot(
      bhist,
      xlim = xlim,
      ylim = ylim,
      col = colors[2],
      border = colors[2],
      add = TRUE
    )
    plot(
      overlap,
      xlim = xlim,
      ylim = ylim,
      col = colors[3],
      border = colors[3],
      add = TRUE
    )
  }

plot_scatterplots <- function(nn, nt, a1, a2, fc, timepointsList, bad1, bad2, x1, x2) {
  #now plot fitted value scatterplots
  
  good <- !bad1 & !bad2
  # TODO: Do these 1, 2 represent hardcoding of # of expected replicates?  If so, need to refactor function
  fit_x1 <- matrix(0, nrow = nn, ncol = nt)
  fit_x2 <- matrix(0, nrow = nn, ncol = nt)
  for (j in 1:nt) {
    fit_x1[, j] <- 2 ^ (a1 + fc * timepointsList[j])
    fit_x2[, j] <- 2 ^ (a2 + fc * timepointsList[j])
  }
  fit_x1[bad1, ] <- 0 #remove missing constructs
  fit_x2[bad2, ] <- 0
  ab1 <- apply(fit_x1, 2, sum)
  ab2 <- apply(fit_x2, 2, sum)
  fit_x1 <- t(t(fit_x1) / ab1) #take ratio to get fitted x_c values
  fit_x2 <- t(t(fit_x2) / ab2)
  fit_x1[!bad1, ] <- log2(fit_x1[!bad1, ])
  fit_x2[!bad2, ] <- log2(fit_x2[!bad2, ])
  fit_x1[bad1, ] <- 0 #remove missing constructs
  fit_x2[bad2, ] <- 0
  
  # TODO: Do these "2:nt"s represent hardcoding of # of expected replicates?  If so, need to refactor function
  #plot experimental values
  for (i in 2:nt) {
    smoothScatter(
      x1[!bad1, c(1, i)],
      xlim = c(-17, -11),
      ylim = c(-17, -11),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      xlab = paste("log2 frequency, d", timepointsList[1], sep = ""),
      ylab = paste("log2 frequency, d", timepointsList[i], sep = ""),
      main = "experimental log2 frequencies, replicate 1"
    )
    abline(0, 1, col = "#000066")
    a <- (x1[!bad1, i] + x1[!bad1, 1]) / 2
    m <- (x1[!bad1, i] - x1[!bad1, 1])
    smoothScatter(
      a,
      m,
      xlim = c(-17, -11),
      ylim = c(-10, 10),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      main = paste(
        "experimental log2 frequencies, replicate 1, d",
        timepointsList[1],
        "-d",
        timepointsList[i],
        sep = ""
      )
    )
    abline(h = 0, col = "#000066")
    yl <- lowess(a, m, f = 0.2)
    lines(yl, col = "red")
  }
  for (i in 2:nt) {
    smoothScatter(
      x2[!bad2, c(1, i)],
      xlim = c(-17, -11),
      ylim = c(-17, -11),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      xlab = paste("log2 frequency, d", timepointsList[1], sep = ""),
      ylab = paste("log2 frequency, d", timepointsList[i], sep = ""),
      main = "experimental log2 frequencies, replicate 2"
    )
    abline(0, 1, col = "#000066")
    a <- (x2[!bad2, i] + x2[!bad2, 1]) / 2
    m <- (x2[!bad2, i] - x2[!bad2, 1])
    smoothScatter(
      a,
      m,
      xlim = c(-17, -11),
      ylim = c(-10, 10),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      main = paste(
        "experimental log2 frequencies, replicate 2, d",
        timepointsList[1],
        "-d",
        timepointsList[i],
        sep = ""
      )
    )
    abline(h = 0, col = "#000066")
    yl <- lowess(a, m, f = 0.2)
    lines(yl, col = "red")
  }
  
  #now plot fitted values
  for (i in 2:nt) {
    smoothScatter(
      fit_x1[!bad1, c(1, i)],
      xlim = c(-17, -11),
      ylim = c(-17, -11),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      xlab = paste("log2 frequency, d", timepointsList[1], sep = ""),
      ylab = paste("log2 frequency, d", timepointsList[i], sep = ""),
      main = "fitted frequencies, replicate 1"
    )
    abline(0, 1, col = "#000066")
    a <- (fit_x1[!bad1, i] + fit_x1[!bad1, 1]) / 2
    m <- (fit_x1[!bad1, i] - fit_x1[!bad1, 1])
    smoothScatter(
      a,
      m,
      xlim = c(-17, -11),
      ylim = c(-10, 10),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      main = paste(
        "fitted log2 frequencies, replicate 1, d",
        timepointsList[1],
        "-d",
        timepointsList[i],
        sep = ""
      )
    )
    abline(h = 0, col = "#000066")
    yl <- lowess(a, m, f = 0.2)
    lines(yl, col = "red")
  }
  for (i in 2:nt) {
    smoothScatter(
      fit_x2[!bad2, c(1, i)],
      xlim = c(-17, -11),
      ylim = c(-17, -11),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      xlab = paste("log2 frequency, d", timepointsList[1], sep = ""),
      ylab = paste("log2 frequency, d", timepointsList[i], sep = ""),
      main = "fitted frequencies, replicate 2"
    )
    abline(0, 1, col = "#000066")
    a <- (fit_x2[!bad2, i] + fit_x2[!bad2, 1]) / 2
    m <- (fit_x2[!bad2, i] - fit_x2[!bad2, 1])
    smoothScatter(
      a,
      m,
      xlim = c(-17, -11),
      ylim = c(-10, 10),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      main = paste(
        "fitted log2 frequencies, replicate 1, d",
        timepointsList[1],
        "-d",
        timepointsList[i],
        sep = ""
      )
    )
    abline(h = 0, col = "#000066")
    yl <- lowess(a, m, f = 0.2)
    lines(yl, col = "red")
  }
  
  #now cross-plots at times 1 2 3..
  for (i in 1:nt) {
    smoothScatter(
      x1[good, i],
      x2[good, i],
      xlim = c(-17, -11),
      ylim = c(-17, -11),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      xlab = paste("replicate 1, d", timepointsList[i], sep = ""),
      ylab = paste("replicate 2, d", timepointsList[i], sep = ""),
      main = paste("comparison of experimental replicates, d", timepointsList[i], sep =
                     "")
    )
    abline(0, 1, col = "#000066")
  }
  for (i in 1:nt) {
    smoothScatter(
      fit_x1[good, i],
      fit_x2[good, i],
      xlim = c(-17, -11),
      ylim = c(-17, -11),
      transformation = function(x) {
        log2(x + 1) ^ .75
      },
      pch = 16,
      cex = 0.4,
      xlab = paste("replicate 1, d", timepointsList[i], sep = ""),
      ylab = paste("replicate 2, d", timepointsList[i], sep = ""),
      main = paste("comparison of fitted replicates, d", timepointsList[i], sep = "")
    )
    abline(0, 1, col = "#000066")
  }
  
}

plotFdrsLeftAndRight <-
  function(pi_mean, pi_null, fdr_left, fdr_right) {
    plot(
      pi_mean,
      fdr_left,
      ylim = c(0, 1),
      xlab = expression(pi["gg'"]),
      ylab = "FDR, left"
    ) #left tail test
    plot(
      pi_mean,
      fdr_right,
      ylim = c(0, 1),
      xlab = expression(pi["gg'"]),
      ylab = "FDR, right"
    ) #right tail test
  }

plotFitToPdf <- function(project, timepointsList, x1, a1, fc, ab1, x2, a2, ab2, nt, pA_pB) {
  pdf(paste(project, ".pdf", sep = ""))
  par(pty = "s", mfrow = c(1, 1))
  # TODO: I don't know what this is doing, but it looks like the method interface has 2 replicates baked right in ...
  plot_fit(timepointsList, x1, a1, fc, ab1, x2, a2, fc, ab2, nt, pA_pB, minfc = 0.10)
  dev.off()
}

plotConstructFitnessesHist <- function(fc, allbad) {
  rge <- range(fc)
  # TODO: where do the 0.005s in the break seq come from?  Do they need to be set in params?
  hh <-
    hist(
      fc[!allbad],
      breaks = seq(rge[1] - 0.005, rge[2] + 0.005, by = 0.005),
      main = "construct fitness",
      col = "grey80",
      border = FALSE,
      xlab = expression(f["c"])
    )
  d <- density(fc[!allbad], bw = 0.01)
  lines(d$x, d$y * sum(hh$counts) * 0.005, col = "black")
  
  rug(fc[!allbad])
  abline(h = 0)
}

plotSingleGeneFitnesses <- function(f, genes) {
  # TODO: since this is always on diagonal, should we refactor to be rug plot?
  rge <- range(f) * 1.1 #with some margin
  plot(
    f[-1],
    f[-1],
    pch = 16,
    cex = 0.6,
    col = "blue",
    xlim = rge,
    ylim = rge,
    xlab = expression(f["g"]),
    ylab = expression(f["g"]),
    main = "single-gene fitness"
  )
  abline(v = 0, h = 0, col = "#000066")
  # TODO: Do I need to worry about these 2s?
  text(
    f[-1],
    f[-1],
    cex = 0.6,
    labels = genes[-1],
    pos = ((rank(f[-1]) %% 2) + 1) * 2,
    offset = 0.2
  )
}

plotRandomVsRealPosteriorProbs <- function(nn, fc, pp_fc, allbad) {
  # nn = number of *good* constructs
  # r = nn random numbers pulled from a uniform distribution with (default) min of 0 and max of 1
  set.seed(1) # TODO: remove after testing.
  r <-
    runif(nn) #TODO: need to seed this to make it deterministic for testing
  fr <-
    fc[r < pp_fc] # errr ... fc for all the constructs where the random ... posterior probability? ... is less than the actual calculated posterior probability for this construct?
  # thus, fr is the distribution of fcs for constructs whose true pp of happening is greater than their random pp of happening--roughly, fcs of constructs whose fc values are more likely than chance
  rge <- range(fc) # range of fcs across all constructs
  
  # TODO: what determines the 0.001 for the break sequence?  Need to set in params?
  plotOverlappingHist(
    fc[!allbad],
    fr,
    breaks = seq(rge[1] - .001, rge[2] + 0.001, by = 0.001),
    xlab = expression(f["c"]),
    ylab = "Frequency"
  )
}