library(dendextend)
library(ComplexHeatmap)
library(beeswarm)
library(survplot)
library(survminer)
library(corrplot)
library(plyr)

# scripts

boop.dir <- function(working.dir, set.dir = FALSE){
  if( !file.exists(working.dir)){
    if( !dir.create(working.dir, recursive = TRUE) ){
      stop(paste("Unable to create directory: ", working.dir, "!\n", sep = ''))
    }
  }
  if(set.dir == TRUE){
    setwd(working.dir)
  }
}

get_split <- function(x, sp = "-", val = 1){
  out.split = unlist(strsplit(x, split = sp, fixed = TRUE))[val]
  return(out.split)
}

split.n.paste <- function(str, sp = ":", start = 1, end = 3, non.sequential = NULL, co = ':'){
  
  boop <- strsplit(x = str, split = sp)
  boop <- unlist(boop)
  if(is.null(non.sequential)){
    boop <- paste(boop[start:end], collapse = co)
  }
  if(!is.null(non.sequential)){
    boop <- paste(boop[non.sequential], collapse = co)
  }
  return(boop)
  
}

classify_based_on_TILs <- function(df){
  df              <- df[!is.na(df$cd8.score.danaher),]
  #df              <- df[!is.na(df$pathology_TILs),]
  df_summary      <- aggregate(df$pathology_TILs, by = list(df$immune_cluster), FUN = median, na.rm=TRUE)
  df_summary$diff <- abs(tmpTILs - df_summary$x)
  return(df_summary$Group.1[which.min(df_summary$diff)])
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      tmp <- cor.test(x = mat[,i], y = mat[,j], ...)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      
      # only "pearson" method provides confidence intervals
      if (!is.null(tmp$conf.int)) {
        lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
        uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
      }
    }
  }
  
  list(
    p = p.mat,
    lowCI = lowCI.mat,
    uppCI = uppCI.mat
  )
}

plot.with.confidence <- function(x, y, x2 = NULL, y2 = NULL, lab = NULL, ci = 0.95, ...){
  
  mini_function <- function(x, y, lab = NULL, ci, ...){
    
    # NA's break it
    z <- c(which(is.na(x)), which(is.na(y)))
    z <- unique(z)
    if(length(z) != 0){
      x <- x[-z]
      y <- y[-z]
      lab <- lab[-z]
    }
    
    df <- data.frame(x, y)
    mod <- lm(y ~ x, data = df)
    
    # get predictions
    new.x <- seq(min(df$x), max(df$x), length.out=length(x))
    #new.y <- seq(min(y), max(y), length.out = length(y))
    preds <- predict(mod, newdata = data.frame(x = new.x), interval = 'confidence', level = ci)
    polygon(c(rev(new.x), new.x), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey90', border = NA)
    
    # line
    abline(mod)
    
    # intervals
    lines(new.x, preds[ ,3], lty = 'dashed', col = 'red')
    lines(new.x, preds[ ,2], lty = 'dashed', col = 'red')
    
    points(x, y, ...)
    
    cor.p <- cor.test(x,y, method = 's')$p.value
    cor.r <- cor.test(x,y, method = 's')$estimate
    
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    ymin <- par("usr")[3]
    ymax <- par("usr")[4]
    
    text(x = 0.7*xmax, y = 0.93*ymax, labels = paste('p: ', formatC(cor.p, format = 'e', digits = 1), '\n', 'rho: ', round(cor.r, digits = 2), sep= ''), pos = 4)
    
    # predict the original values
    real.mod <- predict(mod, level = 0.95, interval = "confidence")
    
    # and are the data points inside the CI or outside of it
    inside.interval <- real.mod[,'lwr'] < y & y < real.mod[,'upr']
    
    if(!is.null(lab)){
      inside.patients  <- lab[inside.interval]
      outside.patients <- lab[!inside.interval] 
    }
    
    if(is.null(lab)){
      inside.patients  <- c()
      outside.patients <- c()
    }
    
    full.model <- list(x, y, lab, new.x, preds, mod, inside.patients, outside.patients)
    names(full.model) <- c('x', 'y', 'lab', 'new.x', 'preds', 'mod', 'inside.patients', 'outside.patients')
    return(full.model)
    
  }
  
  
  # plot empty
  plot(x = c(x, x2), y = c(y, y2), type = 'n', ...)
  
  mini_function(x = x, y = y, lab = lab, ci = ci, ...)
  if(!is.null(x2) & !is.null(y2)){
    mini_function(x = x2, y = y2, lab = lab, ci = ci, ...)
  }
  
}

heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

output.dir           <- '~/Google Drive/SwantonLab/Rimages/test/'
boop.dir(output.dir)

figure.dir           <- paste(output.dir, 'figures-jan15/', sep = '')
boop.dir(figure.dir)




# load some data -- to be saved in an .RData file later #
orig_par <- par()

# ok doing things

all_danaher <- immune.subsets[,c(grep(pattern = 'danah', colnames(immune.subsets), value = TRUE))]
all_danaher <- all_danaher[,c(grep(pattern = 'resid', colnames(all_danaher), value = TRUE, invert = TRUE))]
colnames(all_danaher) <- sapply(colnames(all_danaher), FUN = get_split, sp = '.score', val = 1)

all_danaher_accepted <- all_danaher

# by subtypes
luad_danaher <- all_danaher_accepted[which(substr(rownames(all_danaher_accepted), 1, 8) %in% patient.summary$CRUKid[which(patient.summary$simple.hist %in%  c('LUAD'))]),]
luad_danaher_scaled <- scale(luad_danaher)

lusc_danaher <- all_danaher_accepted[which(substr(rownames(all_danaher_accepted), 1, 8) %in% patient.summary$CRUKid[which(patient.summary$simple.hist == 'LUSC')]),]
lusc_danaher_scaled <- scale(lusc_danaher)

other_danaher <- all_danaher_accepted[which(substr(rownames(all_danaher_accepted), 1, 8) %in% patient.summary$CRUKid[which(patient.summary$simple.hist == 'other')]),]
other_danaher_scaled <- scale(other_danaher)

# empty variables
all_clusters         <- c()
all_classification   <- c()
all_immune_distances <- c()

hmcols         <- colorRampPalette(c('dark blue', 'white','dark red'))(100)

#lusc
hm   <- heatmap.3(t(lusc_danaher_scaled), trace = 'none', hclustfun = function(x) hclust(x, method= 'ward.D2'), dendrogram = 'column', col = hmcols, na.rm=TRUE)
hm_d <- lapply(unique(substr(rownames(lusc_danaher_scaled), 1, 8)), FUN = function(x) {dist(lusc_danaher_scaled[grep(pattern = x, rownames(lusc_danaher_scaled)),,drop=FALSE])} )
names(hm_d)   <- unique(substr(rownames(lusc_danaher_scaled), 1, 8))
all_immune_distances <- c(all_immune_distances, hm_d)
tmp_patients <- unique(substr(rownames(hm$carpet), 1, 8))
tmp_clusters <- cutree(hm$colDendrogram, k = 2)
bp <- boxplot(cd8 ~ tmp_clusters, data = lusc_danaher_scaled, plot=FALSE)
tmp_clusters <- ifelse(tmp_clusters == which.min(bp$stats[3,]), yes = 1, no = 2)
all_clusters <- c(all_clusters, tmp_clusters)
tmp_df       <- as.data.frame(matrix(data = 0, ncol = length(tmp_patients), nrow = length(tmp_clusters), dimnames = list(names(tmp_clusters), tmp_patients)))
for(x in rownames(tmp_df)){
  tmp_df[x,substr(x, 1, 8)] <- tmp_clusters[x]
}
# classify groups
tmp_classify <- apply(tmp_df, MARGIN = 2, table)
tmp_classify_col <- ifelse(lapply(tmp_classify, length) == 3, yes = '#df9b2f', no = ifelse(unlist(lapply(tmp_classify, FUN = function(x) any(grepl(pattern = 1, x = names(x))))), yes = '#30ABDF', no = '#BD2131'))
all_classification <- c(all_classification, tmp_classify_col)
tmp_classify_col <- data.frame(tmp_classify_col, stringsAsFactors = FALSE)
tmp_classify_col$tmp <- sapply(rownames(tmp_classify_col), FUN = function(x){which(colnames(tmp_df) == x)})

lusc_order <- rownames(tmp_classify_col)[with(tmp_classify_col, order(factor(tmp_classify_col, levels = c('#30ABDF', '#df9b2f', '#BD2131')), tmp))]
tmp_df     <- tmp_df[,lusc_order]

# and remove the unnecessary column and re-order the colors
tmp_classify_col$tmp <- NULL
tmp_classify_col <- tmp_classify_col[lusc_order,,drop=FALSE]

# annotate
col_list <- lapply(colnames(tmp_df), FUN = function(x){return(c("1" = '#30ABDF', "0" = "#e9e9e9", "2" = '#BD2131'))})
names(col_list) <- colnames(tmp_df)

ha <- HeatmapAnnotation(tmp_df, show_legend = FALSE, gap = unit(0, 'mm'),show_annotation_name = TRUE, which = 'column', gp = gpar(col = 'white'), annotation_name_gp = gpar(cex = 0.6, col = tmp_classify_col$tmp_classify_col), na_col = 'brown',
                        col = col_list)

hm_lusc <- Heatmap(t(lusc_danaher_scaled), 
                   col = circlize::colorRamp2(c(-3.5, 0, 3.5), c('dark blue', "white", 'dark red')), 
                   clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
                   km = 1,
                   show_column_dend = TRUE, show_row_dend = FALSE, 
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   bottom_annotation = ha,
                   bottom_annotation_height = unit(4, 'in'),
                   row_names_gp = gpar(fontsize = 10))

lusc_cell_type_order <- rev(colnames(lusc_danaher_scaled)[hm$rowInd])


#luad -- all
hm <- heatmap.3(t(luad_danaher_scaled[,lusc_cell_type_order]), trace = 'none', hclustfun = function(x) hclust(x, method= 'ward.D2'), dendrogram = 'column', col = hmcols, na.rm=TRUE, Rowv = FALSE)
#hm <- heatmap.3(t(luad_danaher_scaled), trace = 'none', hclustfun = function(x) hclust(x, method= 'ward.D2'), dendrogram = 'column', col = hmcols, na.rm=TRUE)
hm_d <- lapply(unique(substr(rownames(luad_danaher_scaled), 1, 8)), FUN = function(x) {dist(luad_danaher_scaled[grep(pattern = x, rownames(luad_danaher_scaled)),,drop=FALSE])} )
names(hm_d)   <- unique(substr(rownames(luad_danaher_scaled), 1, 8))
all_immune_distances <- c(all_immune_distances, hm_d)
tmp_patients <- unique(substr(rownames(hm$carpet), 1, 8))
tmp_clusters <- cutree(hm$colDendrogram, k = 2)
bp <- boxplot(cd8 ~ tmp_clusters, data = luad_danaher_scaled, plot=FALSE)
tmp_clusters <- ifelse(tmp_clusters == which.min(bp$stats[3,]), yes = 1, no = 2)
all_clusters <- c(all_clusters, tmp_clusters)
tmp_df       <- as.data.frame(matrix(data = 0, ncol = length(tmp_patients), nrow = length(tmp_clusters), dimnames = list(names(tmp_clusters), tmp_patients)))
for(x in rownames(tmp_df)){
  tmp_df[x,substr(x, 1, 8)] <- tmp_clusters[x]
}
# classify groups
tmp_classify <- apply(tmp_df, MARGIN = 2, table)
tmp_classify_col <- ifelse(lapply(tmp_classify, length) == 3, yes = '#df9b2f', no = ifelse(unlist(lapply(tmp_classify, FUN = function(x) any(grepl(pattern = 1, x = names(x))))), yes = '#30ABDF', no = '#BD2131'))
all_classification <- c(all_classification, tmp_classify_col)
tmp_classify_col <- data.frame(tmp_classify_col, stringsAsFactors = FALSE)
tmp_classify_col$tmp <- sapply(rownames(tmp_classify_col), FUN = function(x){which(colnames(tmp_df) == x)})

luad_order <- rownames(tmp_classify_col)[with(tmp_classify_col, order(factor(tmp_classify_col, levels = c('#30ABDF', '#df9b2f', '#BD2131')), tmp))]
tmp_df     <- tmp_df[,luad_order]

# and remove the unnecessary column and re-order the colors
tmp_classify_col$tmp <- NULL
tmp_classify_col <- tmp_classify_col[luad_order,,drop=FALSE]

# annotate
col_list <- lapply(colnames(tmp_df), FUN = function(x){return(c("1" = '#30ABDF', "0" = "#e9e9e9", "2" = '#BD2131'))})
names(col_list) <- colnames(tmp_df)

ha <- HeatmapAnnotation(tmp_df, show_legend = FALSE, gap = unit(0, 'mm'),show_annotation_name = TRUE, which = 'column', gp = gpar(col = 'white'), annotation_name_gp = gpar(cex = 0.6, col = tmp_classify_col$tmp_classify_col), na_col = 'brown',
                        col = col_list)

hm_luad <- Heatmap(t(luad_danaher_scaled[,lusc_cell_type_order]), 
                   col = circlize::colorRamp2(c(-3.5, 0, 3.5), c('dark blue', "white", 'dark red')), 
                   clustering_method_columns = 'ward.D2', cluster_rows = FALSE, 
                   show_column_dend = TRUE, show_row_dend = FALSE, 
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   bottom_annotation = ha,
                   bottom_annotation_height = unit(4, 'in'),
                   row_names_gp = gpar(fontsize = 10))

#other
hm <- heatmap.3(t(other_danaher_scaled), trace = 'none', hclustfun = function(x) hclust(x, method= 'ward.D2'), dendrogram = 'none', col = hmcols)
hm_d <- lapply(unique(substr(rownames(other_danaher_scaled), 1, 8)), FUN = function(x) {dist(other_danaher_scaled[grep(pattern = x, rownames(other_danaher_scaled)),,drop=FALSE])} )
names(hm_d)   <- unique(substr(rownames(other_danaher_scaled), 1, 8))
all_immune_distances <- c(all_immune_distances, hm_d)
tmp_patients <- unique(substr(rownames(hm$carpet), 1, 8))
tmp_clusters <- cutree(hm$colDendrogram, k = 2)
bp <- boxplot(cd8 ~ tmp_clusters, data = other_danaher_scaled, plot=FALSE)
tmp_clusters <- ifelse(tmp_clusters == which.min(bp$stats[3,]), yes = 1, no = 2)
all_clusters <- c(all_clusters, tmp_clusters)
tmp_df       <- as.data.frame(matrix(data = 0, ncol = length(tmp_patients), nrow = length(tmp_clusters), dimnames = list(names(tmp_clusters), tmp_patients)))
for(x in rownames(tmp_df)){
  tmp_df[x,substr(x, 1, 8)] <- tmp_clusters[x]
}
# classify groups
tmp_classify <- alply(.data = tmp_df, .margins = 2, .fun = table)
names(tmp_classify) <- attr(tmp_classify, "split_labels")$X1
tmp_classify_col <- ifelse(lapply(tmp_classify, length) == 3, yes = '#df9b2f', no = ifelse(unlist(lapply(tmp_classify, FUN = function(x) any(grepl(pattern = 1, x = names(x))))), yes = '#30ABDF', no = '#BD2131'))
all_classification <- c(all_classification, tmp_classify_col)
tmp_classify_col <- data.frame(tmp_classify_col, stringsAsFactors = FALSE)
tmp_classify_col$tmp <- sapply(rownames(tmp_classify_col), FUN = function(x){which(colnames(tmp_df) == x)})

other_order <- rownames(tmp_classify_col)[with(tmp_classify_col, order(factor(tmp_classify_col, levels = c('#30ABDF', '#df9b2f', '#BD2131')), tmp))]
tmp_df     <- tmp_df[,other_order]

# and remove the unnecessary column and re-order the colors
tmp_classify_col$tmp <- NULL
tmp_classify_col <- tmp_classify_col[other_order,,drop=FALSE]

# annotate
col_list <- lapply(colnames(tmp_df), FUN = function(x){return(c("1" = '#30ABDF', "0" = "#e9e9e9", "2" = '#BD2131'))})
names(col_list) <- colnames(tmp_df)

ha <- HeatmapAnnotation(tmp_df, show_legend = FALSE, gap = unit(0, 'mm'),show_annotation_name = TRUE, which = 'column', gp = gpar(col = 'white'), annotation_name_gp = gpar(cex = 0.6, col = tmp_classify_col$tmp_classify_col), na_col = 'brown',
                        col = col_list)

hm_other <- Heatmap(t(other_danaher_scaled), 
                    col = circlize::colorRamp2(c(-3.5, 0, 3.5), c('dark blue', "white", 'dark red')), 
                    clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2', 
                    show_column_dend = TRUE, show_row_dend = FALSE, 
                    show_column_names = FALSE,
                    show_row_names = TRUE,
                    show_heatmap_legend = FALSE,
                    bottom_annotation = ha,
                    bottom_annotation_height = unit(4, 'in'),
                    row_names_gp = gpar(fontsize = 10))

tx100_immune_cluster        <- ifelse(is.na(all_clusters), yes = NA, no = ifelse(all_clusters == 1, yes = 'low', no = 'high'))
names(tx100_immune_cluster) <- names(all_clusters)
tx100_immune_classification <- ifelse(is.na(all_classification), yes = NA, no = ifelse(all_classification == '#BD2131', yes = 'high', no = ifelse(all_classification == '#df9b2f', yes = 'mixed', no = 'low')))
names(tx100_immune_classification) <- names(all_classification)

# add clusters to patient summary
patient.summary$immune_classification <- tx100_immune_classification[match(patient.summary$CRUKid, names(tx100_immune_classification))]
patient.summary.region$immune_cluster <- tx100_immune_cluster[match(patient.summary.region$CRUKidRegion, names(tx100_immune_cluster))]
patient.summary.region$immune_classification <- patient.summary$immune_classification[match(patient.summary.region$CRUKid, patient.summary$CRUKid)]

# salvage the samples without RNAseq...
salvaged_clusters <- c()
for(x in patient.summary.region$CRUKidRegion[is.na(patient.summary.region$cd8.score.danaher)]){
  print(x)
  tmpTILs  <- patient.summary.region[x,'pathology_TILs']
  if(is.na(tmpTILs)){salvaged_clusters <- c(salvaged_clusters, NA)}
  tmp_hist <- patient.summary.region[x,'simple.hist']
  tmp_sub  <- patient.summary.region[x,'subtype']
  if(tmp_hist %in% c('LUSC', 'other', 'LUAD')){
    df     <- patient.summary.region[which(patient.summary.region$simple.hist == tmp_hist),]
    out <- classify_based_on_TILs(df)
  }
  salvaged_clusters <- c(salvaged_clusters, out)
}
names(salvaged_clusters) <- patient.summary.region$CRUKidRegion[is.na(patient.summary.region$cd8.score.danaher)]

salvaged_classification <- sapply(unique(patient.summary.region$CRUKid[is.na(patient.summary.region$cd8.score.danaher)]), FUN = function(x){
  new_clusters <- salvaged_clusters[grep(pattern = x, names(salvaged_clusters))]
  old_clusters <- tx100_immune_cluster[grep(pattern = x, names(tx100_immune_cluster))]
  all_clusters <- c(new_clusters, old_clusters)
  all_clusters <- all_clusters[!is.na(all_clusters)]
  if(length(all_clusters) == 0){return(NA)}
  if(length(table(all_clusters)) == 1){return(names(table(all_clusters)))}
  if(length(table(all_clusters)) > 1){return('mixed')}
})

# combine the two approaches
patient.summary.region$salvaged_clusters       <- salvaged_clusters[match(patient.summary.region$CRUKidRegion, names(salvaged_clusters))]
patient.summary.region$salvaged_classification <- salvaged_classification[match(patient.summary.region$CRUKid, names(salvaged_classification))]
patient.summary$salvaged_classification        <- salvaged_classification[match(patient.summary$CRUKid, names(salvaged_classification))]

# save the original clustering
patient.summary.region$orig_immune_cluster        <- patient.summary.region$immune_cluster
patient.summary.region$orig_immune_classification <- patient.summary.region$immune_classification
patient.summary$orig_immune_classification        <- patient.summary$immune_classification

# and update immune clusters/immune classification
patient.summary.region$immune_cluster[is.na(patient.summary.region$immune_cluster)] <- patient.summary.region$salvaged_clusters[is.na(patient.summary.region$immune_cluster)]
patient.summary$immune_classification <- sapply(1:nrow(patient.summary), FUN = function(x){
  old <- patient.summary$orig_immune_classification[x]
  new <- patient.summary$salvaged_classification[x]
  if(is.na(new)){return(old)}
  else{return(new)}
})
patient.summary.region$immune_classification <- patient.summary$immune_classification[match(patient.summary.region$CRUKid, patient.summary$CRUKid)]

patient.summary$high_immune <- ifelse(patient.summary$immune_classification == 'high', yes = TRUE, no = FALSE)
patient.summary$low_evasion <- ifelse(patient.summary$high_immune == TRUE | (patient.summary$immune_edited == FALSE & patient.summary$HLA_disruption == FALSE), yes = TRUE, no = FALSE)

mutTableRegion$immune_classification <- patient.summary$immune_classification[match(mutTableRegion$CRUKid, patient.summary$CRUKid)]
mutTableRegion$immune_cluster <- patient.summary.region$immune_cluster[match(mutTableRegion$CRUKRegionID, patient.summary.region$CRUKidRegion)]
mutTableAll.together.filtered$immune_classification <- patient.summary$immune_classification[match(mutTableAll.together.filtered$CRUKid, patient.summary$CRUKid)]

keep_to_add$immune_classification <- patient.summary$immune_classification[match(keep_to_add$CRUKid, patient.summary$CRUKid)]
keep_to_add$immune_cluster        <- patient.summary.region$immune_cluster[match(keep_to_add$CRUKRegionID, patient.summary.region$CRUKidRegion)]



# Figure 1A-B
{
  hm_luad
  hm_lusc
  par(orig_par)
}

# Figure 2A
{
  par(mfrow = c(1,1))
  plot(total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'genomic'], total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'immune'], xlab = 'Pairwise genomic distance', ylab = 'Pairwise immune distance', pch = 16, col = total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'color'], ylim = c(0,10), las = 1)
  p.val <- cor.test(total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'genomic'], total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'immune'], method = 's')
  legend(x = 'topright', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\n', 'rho: ', formatC(p.val$estimate, digits = 1, format = 'e'), sep= ''), bty = 'n')
  points(total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'genomic'], total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'immune'], xlab = 'Pairwise genomic distance', ylab = 'Pairwise immune distance', pch = 16, col = total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'color'])
  p.val <- cor.test(total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'genomic'], total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'immune'], method = 's')
  legend(x = 'bottomright', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\n', 'rho: ', formatC(p.val$estimate, digits = 1, format = 'e'), sep= ''), bty = 'n')
}

# Figure 2B-C
{
  par(mfrow = c(1,2))
  
  tmp <- gsub(pattern = 'low', replacement = '#30ABDF', x = gsub(pattern = 'high', replacement = '#BD2131', x = patient.summary.region$immune_cluster[which(patient.summary.region$simple.hist == 'LUAD')]))
  boxplot(shannon.ith ~ factor(immune_classification, levels = c('low', 'mixed', 'high')), data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),], las = 1, ylim = c(0,2.4), main = 'Lung adeno.', ylab = 'Shannon entropy')
  beeswarm(shannon.ith ~ factor(immune_classification, levels = c('low', 'mixed', 'high')), data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),], add = TRUE, pch = 16, pwcol = tmp, cex = 1.6)
  p.val <- wilcox.test(shannon.ith ~ immune_classification == 'low', data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),])$p.value
  legend('bottomright', paste('p -- low-mixed/high : ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n')
  
  tmp <- gsub(pattern = 'low', replacement = '#30ABDF', x = gsub(pattern = 'high', replacement = '#BD2131', x = patient.summary.region$immune_cluster[which(patient.summary.region$simple.hist == 'LUSC')]))
  boxplot(shannon.ith ~ factor(immune_classification, levels = c('low', 'mixed', 'high')), data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUSC'),], las = 1, ylim = c(0,2.4), main = 'Lung squam.', ylab = 'Shannon entropy')
  beeswarm(shannon.ith ~ factor(immune_classification, levels = c('low', 'mixed', 'high')), data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUSC'),], add = TRUE, pch = 16, pwcol = tmp, cex = 1.6)
  p.val <- wilcox.test(shannon.ith ~ immune_classification == 'low', data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUSC'),])$p.value
  legend('bottomright', paste('p -- low-mixed/high : ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n')
}

# Figure 2D
{
  tmp <- patient.summary[,c('ie.all.updated', 'ie.clonal.updated', 'ie.subclonal.updated', 'simple.hist', 'immune_classification')]
  tmp <- tmp[complete.cases(tmp),]
  tmp <- tmp[which(tmp$simple.hist != 'other'),]
  
  par(mfrow = c(1,3))
  
  p.val <- t.test(tmp$ie.clonal.updated[which(tmp$immune_classification == 'high')], tmp$ie.subclonal.updated[which(tmp$immune_classification == 'high')], paired = TRUE)$p.value # less IE in subclones (because clonal IE - subclonal IE < 0)
  stripchart(tmp[which(tmp$immune_classification == 'high'),c('ie.clonal.updated', 'ie.subclonal.updated')], vertical = TRUE, pch= 16, at = c(0.9, 1.2), xlim =c(0.7,1.4), ylim = c(0.7, 1.22), las = 1, group.names = c('Clonal', 'Subclonal'), cex = 1.5, ylab = 'Immunoediting')
  segments(x0 = 0.9, y0 = tmp$ie.clonal.updated[which(tmp$immune_classification == 'high')], x1 = 1.2, y1 = tmp$ie.subclonal.updated[which(tmp$immune_classification == 'high')], col = '#BD2131', lwd = 3)
  legend(legend = paste('paired t-test p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), x = 'topleft', bty = 'n', cex = 0.9)
  
  p.val <- t.test(tmp$ie.clonal.updated[which(tmp$immune_classification == 'mixed')], tmp$ie.subclonal.updated[which(tmp$immune_classification == 'mixed')], paired = TRUE)$p.value # less IE in subclones (because clonal IE - subclonal IE < 0)
  stripchart(tmp[which(tmp$immune_classification == 'mixed'),c('ie.clonal.updated', 'ie.subclonal.updated')], vertical = TRUE, pch= 16, at = c(0.9, 1.2), xlim =c(0.7,1.4), ylim = c(0.7, 1.22), las = 1, group.names = c('Clonal', 'Subclonal'), cex = 1.5, ylab = 'Immunoediting')
  segments(x0 = 0.9, y0 = tmp$ie.clonal.updated[which(tmp$immune_classification == 'mixed')], x1 = 1.2, y1 = tmp$ie.subclonal.updated[which(tmp$immune_classification == 'mixed')], col = '#df9b2f', lwd = 3)
  legend(legend = paste('paired t-test p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), x = 'topleft', bty = 'n', cex = 0.9)
  
  p.val <- t.test(tmp$ie.clonal.updated[which(tmp$immune_classification == 'low')], tmp$ie.subclonal.updated[which(tmp$immune_classification == 'low')], paired = TRUE)$p.value # less IE in subclones (because clonal IE - subclonal IE < 0)
  stripchart(tmp[which(tmp$immune_classification == 'low'),c('ie.clonal.updated', 'ie.subclonal.updated')], vertical = TRUE, pch= 16, at = c(0.9, 1.2), xlim =c(0.7,1.4), ylim = c(0.7, 1.22), las = 1, group.names = c('Clonal', 'Subclonal'), cex = 1.5, ylab = 'Immunoediting')
  segments(x0 = 0.9, y0 = tmp$ie.clonal.updated[which(tmp$immune_classification == 'low')], x1 = 1.2, y1 = tmp$ie.subclonal.updated[which(tmp$immune_classification == 'low')], col = '#30ABDF', lwd = 3)
  legend(legend = paste('paired t-test p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), x = 'topleft', bty = 'n', cex = 0.9)
}

# Figure 2F
{
  tmp <- mutTableRegion[which(mutTableRegion$mutation_id %in% mutTableAll.together.filtered$mutation_id[!is.na(mutTableAll.together.filtered$CPNremove)]),]
  tmp <- tmp[which(tmp$nonsynonymous == TRUE),]
  tmp <- tmp[which(tmp$DriverMut == FALSE),]
  
  tmp <- tmp[,c('CRUKRegionID', 'mutation_id', 'count', 'strong.count', 'any.HLA.loss', 'immune_classification', 'immune_cluster', 'nonsynonymous', 'DriverMut', 'CPNremove')]
  to_add <- keep_to_add[,c('CRUKRegionID', 'mutation_id', 'count', 'strong.count', 'any.HLA.loss', 'immune_classification', 'immune_cluster', 'nonsynonymous', 'DriverMut', 'CPNremove')]
  
  tmp <- rbind(tmp, to_add)

  par(mfrow = c(2,1))

  tmp_nodup  <- tmp[!duplicated(tmp$mutation_id),]
  to_barplot <- sort(table(substr(tmp_nodup$CRUKRegionID, 1, 8)[which(tmp_nodup$count == 1 & tmp_nodup$CPNremove == TRUE)]), decreasing=TRUE)
  barplot(to_barplot, las = 1, names = NA, ylim = c(0,40), col = 'black', ylab = 'Number CN loss neo', xlab = 'Tumor sample')
  
  props <- to_barplot/(patient.summary$ClonalNeo[match(names(to_barplot), patient.summary$CRUKid)])
  names(props) <- NA
  barplot(rbind(props, 0.5-props), las = 1, ylim = c(0,0.5), col = c('black', 'lightgrey'), ylab = '% CN loss neo', xlab = 'Tumor sample', border = NA)
  
}

# Figure 2G
{

  tmp <- mutTableRegion[which(mutTableRegion$mutation_id %in% mutTableAll.together.filtered$mutation_id[!is.na(mutTableAll.together.filtered$CPNremove)]),]
  tmp <- tmp[which(tmp$nonsynonymous == TRUE),]
  tmp <- tmp[which(tmp$DriverMut == FALSE),]
  
  tmp <- tmp[,c('regionMutID', 'mutation_id', 'count', 'strong.count', 'any.HLA.loss', 'immune_classification', 'immune_cluster', 'nonsynonymous', 'DriverMut', 'CPNremove')]
  to_add <- keep_to_add[,c('regionMutID', 'mutation_id', 'count', 'strong.count', 'any.HLA.loss', 'immune_classification', 'immune_cluster', 'nonsynonymous', 'DriverMut', 'CPNremove')]
  
  tmp <- rbind(tmp, to_add)
  
  {par(mfrow = c(1,1))
    ests   <- c()
    p.vals <- c()
    ors1   <- c()
    ors2   <- c()
    neo_ratio <- c()
    non_neo_ratio <- c()  
    
    tmp2 <- table(tmp$count, tmp$CPNremove)
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    neo_ratio <- c(neo_ratio, paste(as.character(tmp2['1','TRUE']), '/', as.character(tmp2['1', 'FALSE'] + tmp2['1', 'TRUE']), sep = ''))
    non_neo_ratio <- c(non_neo_ratio, paste(as.character(tmp2['0','TRUE']), '/', as.character(tmp2['0', 'FALSE'] + tmp2['0', 'TRUE']), sep = ''))
    
    tmp2 <- table(tmp$count[which(tmp$immune_cluster %in% c('high'))], tmp$CPNremove[which(tmp$immune_cluster %in% c('high'))])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    neo_ratio <- c(neo_ratio, paste(as.character(tmp2['1','TRUE']), '/', as.character(tmp2['1', 'FALSE'] + tmp2['1', 'TRUE']), sep = ''))
    non_neo_ratio <- c(non_neo_ratio, paste(as.character(tmp2['0','TRUE']), '/', as.character(tmp2['0', 'FALSE'] + tmp2['0', 'TRUE']), sep = ''))
    
    tmp2 <- table(tmp$count[which(tmp$immune_cluster %in% c('low'))], tmp$CPNremove[which(tmp$immune_cluster %in% c('low'))])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    neo_ratio <- c(neo_ratio, paste(as.character(tmp2['1','TRUE']), '/', as.character(tmp2['1', 'FALSE'] + tmp2['1', 'TRUE']), sep = ''))
    non_neo_ratio <- c(non_neo_ratio, paste(as.character(tmp2['0','TRUE']), '/', as.character(tmp2['0', 'FALSE'] + tmp2['0', 'TRUE']), sep = ''))
    
    ests <- ests
    ors1 <- ors1
    ors2 <- ors2
    
    bp <- barplot(ests, col = c('#00441b', '#BD2131', '#30ABDF'), border = FALSE, las = 1, cex.axis = 0.9, cex.names = 0.7, ylab = 'OR -- neo in region of CPN loss', ylim = c(0,3), names = c('All', 'High', 'Low'), axes = FALSE)
    axis(side = 2, at = c(0,1,2,3), labels = c(0,1,2,3), las = 2)
    mtext(text = paste('p = ', formatC(p.vals, digits = 1, format = 'e'), sep = ''), side = 1, line = 2, cex = 0.7, at = bp[,1])
    segments(x0 = bp[1,1], x1 = bp[1,1], y0 = ors1[1], y1 = ors2[1])
    segments(x0 = bp[2,1], x1 = bp[2,1], y0 = ors1[2], y1 = ors2[2])
    segments(x0 = bp[3,1], x1 = bp[3,1], y0 = ors1[3], y1 = ors2[3])
    abline(h = 1, lty = 2, col = 'black')
    
    mtext(text = neo_ratio, side = 1, line = 3, cex = 0.7, at = bp[,1])
    mtext(text = non_neo_ratio, side = 1, line = 4, cex = 0.7, at = bp[,1])
  }
}

# Figure 2H
{
  tmp <- patient.summary[,c('ie.all.updated', 'ie.clonal.updated', 'ie.subclonal.updated', 'simple.hist', 'immune_classification', 'anyCPN_remove')]
  tmp <- tmp[complete.cases(tmp),]
  tmp <- tmp[which(tmp$simple.hist != 'other'),]

  par(mfrow = c(1,2))
  
  p.val <- t.test(tmp$ie.clonal.updated[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==TRUE)], tmp$ie.subclonal.updated[which(tmp$immune_classification == 'low' & tmp$anyCPN_remove == TRUE)], paired = TRUE)$p.value # less IE in subclones (because clonal IE - subclonal IE < 0)
  stripchart(tmp[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==TRUE),c('ie.clonal.updated', 'ie.subclonal.updated')], vertical = TRUE, pch= 16, at = c(0.9, 1.2), xlim =c(0.7,1.4), ylim = c(0.7, 1.22), las = 1, group.names = c('Clonal', 'Subclonal'), cex = 1.5, ylab = 'Immunoediting')
  segments(x0 = 0.9, y0 = tmp$ie.clonal.updated[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==TRUE)], x1 = 1.2, y1 = tmp$ie.subclonal.updated[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==TRUE)], col = '#30ABDF', lwd = 3)
  legend(legend = paste('paired t-test p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), x = 'topleft', bty = 'n', cex = 0.9)
  
  p.val <- t.test(tmp$ie.clonal.updated[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==FALSE)], tmp$ie.subclonal.updated[which(tmp$immune_classification == 'low' & tmp$anyCPN_remove == FALSE)], paired = TRUE)$p.value # less IE in subclones (because clonal IE - subclonal IE < 0)
  stripchart(tmp[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==FALSE),c('ie.clonal.updated', 'ie.subclonal.updated')], vertical = TRUE, pch= 16, at = c(0.9, 1.2), xlim =c(0.7,1.4), ylim = c(0.7, 1.22), las = 1, group.names = c('Clonal', 'Subclonal'), cex = 1.5, ylab = 'Immunoediting')
  segments(x0 = 0.9, y0 = tmp$ie.clonal.updated[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==FALSE)], x1 = 1.2, y1 = tmp$ie.subclonal.updated[which(tmp$immune_classification == 'low'& tmp$anyCPN_remove==FALSE)], col = '#30ABDF', lwd = 3)
  legend(legend = paste('paired t-test p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), x = 'topleft', bty = 'n', cex = 0.9)
}

# Figure 3A
{
  layout(mat = matrix(c(1,2,3), nrow = 3), heights = c(0.7, 0.2, 0.1))
  par(mar = c(2,4,2,2))
  tmp <- patient.summary[,c('ClonalNeo', 'SubclonalNeo', 'TotalExpressedNeo', 'ClonalExpressedNeo', 'SubclonalExpressedNeo', 'immune_classification', 'simple.hist')]
  tmp$simple.hist <- factor(tmp$simple.hist, levels = c('LUAD', 'LUSC', 'other'))
  tmp <- tmp[which(tmp$simple.hist != 'other'),]
  tmp <- tmp[!is.na(tmp$TotalExpressedNeo),]
  tmp <- tmp[-which(tmp$ClonalNeo + tmp$SubclonalNeo == 0),]
  tmp <- tmp[order(tmp$ClonalNeo, decreasing = TRUE),]
  tmp <- tmp[order(tmp$simple.hist),]
  tmp$ClonalExpressedNeoDif <- tmp$ClonalNeo-tmp$ClonalExpressedNeo
  tmp$SubclonalExpressedNeoDif <- tmp$SubclonalNeo-tmp$SubclonalExpressedNeo
  ord <- rownames(tmp)
  
  bp <- barplot(unname(t(tmp[,c('ClonalExpressedNeo', 'ClonalExpressedNeoDif')])), width = 0.5, space = 2, col = c(c('#3182bd', '#3182bd90')), border = FALSE, las = 2, cex.axis = 0.9, ylab = 'Number Neo', ylim = c(0,800))
  barplot(t(tmp[,c('SubclonalExpressedNeo', 'SubclonalExpressedNeoDif')]), width = 0.5, space = c(3,rep(2, nrow(tmp)-1)), col = c(c('#de2d26', '#de2d2690')), border = FALSE, add = TRUE, axes = FALSE, axisnames = FALSE)
  legend('topright', legend = c('Clonal Neo', 'Clonal Neo (exp)', 'Subclonal Neo', 'Subclonal Neo (exp)'), fill = c('#3182bd90', '#3182bd', '#de2d2690', '#de2d26'), border = FALSE, bty = 'n', ncol = 2, cex = 0.7)
  
  tmp <- mutTableAll.together.filtered[which(mutTableAll.together.filtered$count == 1),]
  tmp <- tmp[!is.na(tmp$ITHState_RNA),]
  tmp <- tmp[which(tmp$PyCloneClonal == 'C'),]
  to_plot <- as.data.frame.matrix(table(tmp$CRUKid, tmp$ITHState_RNA), stringsAsFactors = FALSE)
  to_plot <- to_plot/(rowSums(to_plot))
  to_plot <- to_plot[ord,]
  
  barplot(rbind(to_plot[,1], 1-to_plot[,1]), col = c('black', 'lightgrey'), border = FALSE, width = 0.5, axisnames = FALSE, las = 2, cex.names = 0.7, ylab = 'Fraction of clonal neoantigens ubiquitously expressed')
  tmp_col <- ifelse(patient.summary$simple.hist[match(rownames(to_plot), patient.summary$CRUKid)] == 'LUAD', yes = 'dark blue', no = ifelse(patient.summary$simple.hist[match(rownames(to_plot), patient.summary$CRUKid)] == 'LUSC', yes = 'dark red', no = 'grey'))
  barplot(height = rep(1, nrow(to_plot)), width = 0.5, axisnames = FALSE, col = tmp_col, border = FALSE, axes = FALSE)

  to_plot$immune_classification <- patient.summary$immune_classification[match(rownames(to_plot), patient.summary$CRUKid)] 
  
  par(orig_par)
}

# Figure 3B
{
  
  tmp <- mutTableAll.together.filtered[which(mutTableAll.together.filtered$nonsynonymous == FALSE),]
  tmp <- tmp[!is.na(tmp$ITHState_RNA),]
  tmp <- tmp[which(tmp$PyCloneClonal == 'C'),]
  tmp <- tmp[which(tmp$DriverMut == FALSE),]
  not_to_plot <- as.data.frame.matrix(table(tmp$CRUKid, tmp$ITHState_RNA), stringsAsFactors = FALSE)
  not_to_plot <- not_to_plot/(rowSums(not_to_plot))
  not_to_plot <- not_to_plot[ord,]
  
  patient.summary$fraction_clonal_neo_ubiquitously_expressed <- to_plot[,1][match(patient.summary$CRUKid, rownames(to_plot))]
  patient.summary$immune_classification <- factor(patient.summary$immune_classification, levels = c('low', 'mixed', 'high'))
  par(mfrow = c(1,1))
  boxplot(fraction_clonal_neo_ubiquitously_expressed ~ immune_classification, data = patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),], ylab = 'Fraction of clonal neoantigens UBIQUITOUSLY expressed', las = 1, ylim = c(0,1), outpch = NA)
  beeswarm(fraction_clonal_neo_ubiquitously_expressed ~ immune_classification, data = patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),], method = 'swarm', corral = 'wrap', add = TRUE, pch = 16, pwcol = ifelse(patient.summary$simple.hist[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC'))] == 'LUAD', yes = 'darkblue', no = ifelse(patient.summary$simple.hist[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC'))] == 'LUSC', yes = 'firebrick4', no = 'lightgrey')))
  p.val <- wilcox.test(fraction_clonal_neo_ubiquitously_expressed ~ immune_classification == 'low', data = patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),])$p.value
  legend(x = 'topright', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
  
}

# Figure 3C
{
  
  par(mfrow = c(1,1))
  tmp <- mutTableAll.together.filtered[which(mutTableAll.together.filtered$ITHState_RNA %in% c(1,2,3, 'FALSE')),] 

  {
    par(mfrow = c(1,1))
    ests   <- c()
    p.vals <- c()
    ors1    <- c()
    ors2    <- c()
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == TRUE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == FALSE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == TRUE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == FALSE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == TRUE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == FALSE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == TRUE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == FALSE),'count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    ests <- 1/ests
    ors1 <- 1/ors1
    ors2 <- 1/ors2
    
    bp <- barplot(ests, col = c(sapply(X = c('#00441b', '#BD2131', '#df9b2f', '#30ABDF'), FUN = function(x) rep(x, 2))), border = FALSE, las = 1, cex.axis = 0.9, cex.names = 0.7, ylab = 'OR -- neoantigen expression lost in RNA', ylim = c(0,1.5), names = '', axes = FALSE, space = rep(c(0.6,0.2),4))
    axis(side = 2, at = c(0,0.5, 1.0, 1.5), labels = c(0, 0.5, 1.0, 1.5), las = 1, cex.axis = 0.8)
    segments(x0 = bp[1,1], x1 = bp[1,1], y0 = ors1[1], y1 = ors2[1])
    segments(x0 = bp[2,1], x1 = bp[2,1], y0 = ors1[2], y1 = ors2[2])
    segments(x0 = bp[3,1], x1 = bp[3,1], y0 = ors1[3], y1 = ors2[3])
    segments(x0 = bp[4,1], x1 = bp[4,1], y0 = ors1[4], y1 = ors2[4])
    segments(x0 = bp[5,1], x1 = bp[5,1], y0 = ors1[5], y1 = ors2[5])
    segments(x0 = bp[6,1], x1 = bp[6,1], y0 = ors1[6], y1 = ors2[6])
    segments(x0 = bp[7,1], x1 = bp[7,1], y0 = ors1[7], y1 = ors2[7])
    segments(x0 = bp[8,1], x1 = bp[8,1], y0 = ors1[8], y1 = ors2[8])
    
    mtext(text = paste('p = ', round(p.vals, digits = 3), sep = ''), side = 1, line = 2, cex = 0.7, at = bp[,1])
    mtext(text = c('+HLA LOH', '-HLA LOH'), side = 1, line = 1, cex = 0.7, at = bp[,1])
    tmp_names <- factor(rep(c('All', 'High', 'Hetero.', 'Low'), each = 2), levels = c('All', 'High', 'Hetero.', 'Low'))
    mtext(text = c('All', 'High', 'Hetero.', 'Low'), side = 1, line = 3, cex = 0.7, at = aggregate(bp, by = list(tmp_names), mean)$V1)
    abline(h = 1, lty = 2, col = 'black')
    
  }
  
}

# Figure 3D
{

  {par(mfrow = c(1,1))
    ests   <- c()
    p.vals <- c()
    ors1    <- c()
    ors2    <- c()
    
    tmp2 <- table(mutTableAll.together.filtered$tcga_gene_minimally_expressed_95[which(mutTableAll.together.filtered$nonsynonymous == TRUE)], mutTableAll.together.filtered$count[which(mutTableAll.together.filtered$nonsynonymous == TRUE)])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(mutTableAll.together.filtered$tcga_gene_minimally_expressed_95[which(mutTableAll.together.filtered$nonsynonymous == TRUE & mutTableAll.together.filtered$immune_classification == 'high')], mutTableAll.together.filtered$count[which(mutTableAll.together.filtered$nonsynonymous == TRUE & mutTableAll.together.filtered$immune_classification == 'high')])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(mutTableAll.together.filtered$tcga_gene_minimally_expressed_95[which(mutTableAll.together.filtered$nonsynonymous == TRUE & mutTableAll.together.filtered$immune_classification == 'mixed')], mutTableAll.together.filtered$count[which(mutTableAll.together.filtered$nonsynonymous == TRUE & mutTableAll.together.filtered$immune_classification == 'mixed')])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(mutTableAll.together.filtered$tcga_gene_minimally_expressed_95[which(mutTableAll.together.filtered$nonsynonymous == TRUE & mutTableAll.together.filtered$immune_classification == 'low')], mutTableAll.together.filtered$count[which(mutTableAll.together.filtered$nonsynonymous == TRUE & mutTableAll.together.filtered$immune_classification == 'low')])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    bp <- barplot(ests, col = c('#00441b', '#BD2131', '#df9b2f', '#30ABDF'), border = FALSE, las = 1, cex.axis = 0.9, cex.names = 0.7, ylab = 'OR -- neoantigens depleted in (TCGA) expressed genes', ylim = c(0,1), names = c('All', 'High', 'Hetero.', 'Low'), axes = FALSE)
    axis(side = 2, at = c(0,0.5, 1.0), labels = c(0, 0.5, 1.0), las = 1, cex.axis = 0.8)
    segments(x0 = bp[1,1], x1 = bp[1,1], y0 = ors1[1], y1 = ors2[1])
    segments(x0 = bp[2,1], x1 = bp[2,1], y0 = ors1[2], y1 = ors2[2])
    segments(x0 = bp[3,1], x1 = bp[3,1], y0 = ors1[3], y1 = ors2[3])
    segments(x0 = bp[4,1], x1 = bp[4,1], y0 = ors1[4], y1 = ors2[4])
    mtext(text = paste('p = ', formatC(p.vals, digits = 1, format = 'e'), sep = ''), side = 1, line = 2, cex = 0.7, at = bp[,1])
    abline(h = 1, lty = 2, col = 'black')
  }
}

# Figure 4A-B
{
  for(cancer in c('LUAD', 'LUSC')){
    if(cancer == 'LUAD'){
      colOrder <- c("CRUK0020", "CRUK0016", "CRUK0039", "CRUK0051", "CRUK0027", "CRUK0060", "CRUK0024", "CRUK0035", "CRUK0047", "CRUK0052", "CRUK0026", "CRUK0008", "CRUK0038", "CRUK0045", "CRUK0031", "CRUK0034", "CRUK0006", "CRUK0033", "CRUK0030", "CRUK0015", "CRUK0014", "CRUK0019", "CRUK0007", "CRUK0040", "CRUK0017", "CRUK0001", "CRUK0009", "CRUK0003", "CRUK0029", "CRUK0061", "CRUK0044", "CRUK0004", "CRUK0012", "CRUK0005", "CRUK0041", "CRUK0018", "CRUK0046", "CRUK0048", "CRUK0032", "CRUK0010", "CRUK0013", "CRUK0028", "CRUK0002", "CRUK0022", "CRUK0037", "CRUK0055", "CRUK0023", "CRUK0036", "CRUK0042", "CRUK0053", "CRUK0057", "CRUK0025", "CRUK0050", "CRUK0021", "CRUK0058", "CRUK0011", "CRUK0043", "CRUK0049", "CRUK0054", "CRUK0056", "CRUK0059")
    }
    if(cancer == 'LUSC'){
      colOrder <- c(c("CRUK0086","CRUK0074","CRUK0068","CRUK0079","CRUK0069","CRUK0064","CRUK0083","CRUK0072","CRUK0073","CRUK0081","CRUK0065","CRUK0078","CRUK0067","CRUK0070","CRUK0062","CRUK0082","CRUK0076","CRUK0071","CRUK0066","CRUK0084","CRUK0090","CRUK0075","CRUK0063","CRUK0089","CRUK0085","CRUK0088","CRUK0087","CRUK0092","CRUK0077","CRUK0080","CRUK0093","CRUK0091"))
    }
    
    tmp_df <- data.frame(colnames(all_condensed_alterations)[which(substr(colnames(all_condensed_alterations), 1, 8) %in% colOrder)], stringsAsFactors = FALSE)
    tmp_df$SampleID <- substr(tmp_df[,1], 1, 8)
    tmp_df$SampleID <- factor(tmp_df$SampleID, levels = colOrder)
    tmp_df <- tmp_df[order(tmp_df$SampleID),]
    
    condensed_alterations <- all_condensed_alterations[c('HLALOH', 'APC'),tmp_df[,1]]
    
    
    # add in some clinical things
    tmp_clinical <- matrix(NA, nrow=7, ncol=ncol(condensed_alterations), dimnames = list(c('ClonalNeo', 'SubclonalNeo', 'immune', 'PFS', 'packyears', 'IE', 'low_evasion'), colnames(condensed_alterations)))
    tmp_clinical['ClonalNeo',]     <- patient.summary.region$ClonalNeo[match(colnames(tmp_clinical), patient.summary.region$CRUKidRegion)]
    tmp_clinical['SubclonalNeo',]     <- patient.summary.region$SubclonalNeo[match(colnames(tmp_clinical), patient.summary.region$CRUKidRegion)]
    tmp_clinical['immune',]    <- patient.summary.region$immune_cluster[match(colnames(tmp_clinical), patient.summary.region$CRUKidRegion)]
    tmp_clinical['IE',]        <- patient.summary$immune_edited[match(substr(colnames(tmp_clinical), 1, 8), patient.summary$CRUKid)]
    tmp_clinical['PFS',]          <- patient.summary$allan_final_dsf_time[match(substr(colnames(tmp_clinical), 1, 8), patient.summary$CRUKid)]
    cens                          <- patient.summary$allan_final_dsf[match(substr(colnames(tmp_clinical), 1, 8), patient.summary$CRUKid)]
    tmp_clinical['packyears',]    <- patient.summary$packyears[match(substr(colnames(tmp_clinical), 1, 8), patient.summary$CRUKid)]
    tmp_clinical['low_evasion',]  <- patient.summary$low_evasion[match(substr(colnames(tmp_clinical), 1, 8), patient.summary$CRUKid)]
    
    color.gradient <- function(x, colors=c("dark blue","white","dark red"), colsteps=81) {
      return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
    }
    
    
    pfs_times <- seq(from = 34, to = 1275, by = 1)
    pfs_times_col_1  <- color.gradient(pfs_times, colors = c('#980043', '#d4b9da'), colsteps = 1242)
    pfs_times_col_0  <- color.gradient(pfs_times, colors = c('#ccece6', '#006d2c'), colsteps = 1242)
    
    # packyears
    tmp_smoke        <- as.numeric(sort(unique(patient.summary[,'packyears'])))*100
    tmp_smoke        <- tmp_smoke[!is.na(tmp_smoke)]
    packyearcols <- color.gradient(tmp_smoke, colors = c('lightblue', 'black'))
    tmp_smoke        <- tmp_smoke/100
    
    # plot
    layout(matrix(1:(nrow(condensed_alterations) + nrow(tmp_clinical)), ncol = 1))
    par(mar = c(2,2,1,2), oma = c(4,4,4,4))
    tmp_acceptedRegions_df <- acceptedRegions_df[match(colnames(condensed_alterations), acceptedRegions_df$CRUKidRegion),]
    
    rbPal_nsMut     <- colorRampPalette(c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))
    
    tmp_clinical_colors <- tmp_clinical
    tmp_clinical_colors['ClonalNeo',]    <- rbPal_nsMut(10)[as.numeric(cut(as.numeric(tmp_clinical['ClonalNeo',]),breaks=10))]
    tmp_clinical_colors['SubclonalNeo',] <- rbPal_nsMut(10)[as.numeric(cut(as.numeric(tmp_clinical['SubclonalNeo',]),breaks=10))]
    tmp_clinical_colors['immune',] <- ifelse(is.na(tmp_clinical_colors['immune',]), yes = 'darkgrey', no = ifelse(tmp_clinical_colors['immune',] == 'high', yes = '#BD2131', no = '#30ABDF'))
    tmp_clinical_colors['PFS',which(cens == 0)] <- pfs_times_col_0[match(tmp_clinical['PFS',which(cens == 0)], pfs_times)]
    tmp_clinical_colors['PFS',which(cens == 1)] <- pfs_times_col_1[match(tmp_clinical['PFS',which(cens == 1)], pfs_times)]
    tmp_clinical_colors['packyears',]            <- packyearcols[match(as.numeric(tmp_clinical_colors['packyears',]), tmp_smoke)]
    tmp_clinical_colors['IE',] <- ifelse(tmp_clinical['IE',] == 'TRUE', yes = '#6a3d9a', no = '#cab2d6')
    tmp_clinical_colors['low_evasion',] <- ifelse(tmp_clinical['low_evasion',] == 'TRUE', yes = '#2356A7', no = '#BE1E2D')
    tmp_clinical_colors['low_evasion',] <- ifelse(is.na(tmp_clinical_colors['low_evasion',]), yes = 'darkgrey', no = tmp_clinical_colors['low_evasion',])
    
    bp <-  barplot(tmp_clinical[1:2,],
                   border = NA,axes = FALSE, axisnames = FALSE,
                   col = c('#377eb8', '#e41a1c'), ylim = c(0,1100),
                   space = unlist(sapply(tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)], FUN = function(x){c(2,rep(0.2,(x-1)))}))
    )
    axis(side = 2, at = c(0,500,1000), labels = c(0,500, 1000), line = 1, las = 1)
    mtext(text = 'Neo', side=2,line = -25, cex = 5, las = 2, outer = FALSE)
    
    sapply(3:7, FUN = function(x){
      bp <-  barplot(rep(1,ncol(tmp_clinical_colors)),
                     border = NA,axes = FALSE,
                     col = tmp_clinical_colors[x,], 
                     space = unlist(sapply(tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)], FUN = function(x){c(2,rep(0.2,(x-1)))}))
      )
    })
    
    sapply(1:nrow(condensed_alterations), FUN = function(x){
      bp <-  barplot(rep(1,ncol(condensed_alterations)),border = NA,axes = FALSE,
                     col = ifelse(as.numeric(condensed_alterations[x,]) == 1, yes = '#008000', no = ifelse(as.numeric(condensed_alterations[x,]) == -1, yes = 'blue', no = 'lightgrey')), 
                     space = unlist(sapply(tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)], FUN = function(x){c(2,rep(0.2,(x-1)))}))
      )
      mtext(text = rownames(condensed_alterations)[x], side=2,line = -25, cex = 5, las = 2, outer = FALSE, col = 'black')
    })
    bp <-  barplot(rep(1,ncol(condensed_alterations)),border = NA,axes = FALSE,space = unlist(sapply(tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)], FUN = function(x){c(2,rep(0.2,(x-1)))})), plot = FALSE)
    nr          <- tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)]
    names(nr)   <- unique(tmp_acceptedRegions_df$CRUKid)
    text_start  <- sapply(names(nr), FUN = function(x) {bp[(sum(nr[1:which(names(nr) == x)]) - nr[which(names(nr) == x)] + 1),1]} )
    mtext(text = patient.summary$CRUKid[match(substr(names(nr), 3, nchar(names(nr))), patient.summary$CRUKid)], side=1,line = -1, cex = 5, las = 2, outer = FALSE, at = text_start, col = 'black')
  }
  
  par(orig_par)
}

# Figure 4C
{
  tmp <- patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),]

  tmp_low_clonal <- tmp[which(tmp$high_clonal == FALSE),]
  
  par(mfrow = c(1,3))
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ low_evasion, 
           data = tmp, 
           main = "LUAD&LUSC\nlow immune evasion", 
           xlab = 'Time (days)', 
           ylab = "PFS",
           col = c('darkred', 'darkblue'),
           lwd = 1, las = 1, cex.axis = 0.8)
  
  
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ low_evasion, 
           data = tmp_low_clonal, 
           main = "LUAD&LUSC\nlow immune evasion, low clonal neo", 
           xlab = 'Time (days)', 
           ylab = "PFS", 
           col = c('darkred', 'darkblue'),
           lwd = 1, las = 1, cex.axis = 0.8)
  
  tmp$low_evasion_clonal <- ifelse(tmp$low_evasion == TRUE | tmp$high_clonal == TRUE, yes =TRUE, no = FALSE)
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ low_evasion_clonal, 
           data = tmp, 
           main = "LUAD&LUSC\nlow immune evasion, high clonal neo", 
           xlab = 'Time (days)', 
           ylab = "PFS", 
           col = c('darkred', 'darkblue'),
           lwd = 1, las = 1, cex.axis = 0.8)
}

# Figure S3
{
  par(mfrow = c(1,1))
  boxplot(cd8.score.danaher ~ as.character(simple.hist), data = patient.summary.region[which(patient.summary.region$simple.hist %in% c('LUAD', 'LUSC')),], las = 1, ylab = 'Danaher CD8+ Score')
  beeswarm(cd8.score.danaher ~ as.character(simple.hist), data = patient.summary.region[which(patient.summary.region$simple.hist %in% c('LUAD', 'LUSC')),], add = TRUE, pch = 16)
  p.val <- wilcox.test(cd8.score.danaher ~ simple.hist, patient.summary.region[which(patient.summary.region$simple.hist %in% c('LUAD', 'LUSC')),])$p.value
  legend(x = 'topright', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
}



# Figure S5A
{
  par(mfrow = c(1,1))
  tmp <- patient.summary.region
  boxplot(pathology_TILs ~ factor(orig_immune_cluster, levels = c('low', 'high')), data = tmp, las = 1, ylab = 'pathology estimated TILs', xlab = 'immune cluster')
  beeswarm(pathology_TILs ~ factor(orig_immune_cluster, levels = c('low', 'high')), data = tmp, pch = 16, add = TRUE, cex = 1.5, pwcol = ifelse(is.na(tmp$orig_immune_cluster), yes = 'light grey', no = ifelse(tmp$orig_immune_cluster == 'high', yes = '#BD2131', no = '#30ABDF')))
  p.val <- wilcox.test(pathology_TILs ~ orig_immune_cluster, data = tmp)$p.value
  legend(x = 'topright', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
}

# Figure S5B
{
  par(mfrow = c(1,1))
  tmp_all_danaher_accepted <- all_danaher_accepted[,c('cd8', 'cd4', 'bcell', 'cd45', 'cyto' ,'dend', 'mast', 'nkcd56dim', 'nk', 'tcells', 'th1')]
  tmp_all_danaher_accepted$MDSC_jiang   <- patient.summary.region$MDSC_jiang[match(rownames(tmp_all_danaher_accepted), patient.summary.region$CRUKidRegion)]
  tmp_all_danaher_accepted$CAF_jiang    <- patient.summary.region$CAF_jiang[match(rownames(tmp_all_danaher_accepted), patient.summary.region$CRUKidRegion)]
  tmp_all_danaher_accepted$TAM.M2_jiang <- patient.summary.region$TAM.M2_jiang[match(rownames(tmp_all_danaher_accepted), patient.summary.region$CRUKidRegion)]
  corr <- cor.mtest(tmp_all_danaher_accepted, method = 's')
  M    <- cor(tmp_all_danaher_accepted, method = 's')
  corrplot(M, method = 'square', p.mat = corr$p, diag = TRUE, lowCI.mat = corr$lowCI, uppCI.mat = corr$uppCI, type = 'upper')
}

# Figure S5C
{
  par(mfrow = c(1,2))
  boxplot(TAM.M2_jiang ~ immune_cluster, data = patient.summary.region, las = 1, ylab = 'TAM M2 (Jiang)', xlab = 'Immune classification', out.pch = NA)
  beeswarm(TAM.M2_jiang ~ immune_cluster, data = patient.summary.region, pch = 16, add = TRUE, col = 'black', method = 'swarm', corral = 'wrap')
  p.val <- wilcox.test(TAM.M2_jiang ~ immune_cluster, data = patient.summary.region)$p.value
  legend(x = 'bottomright', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
  
  boxplot(MDSC_jiang ~ immune_cluster, data = patient.summary.region, las = 1, ylab = 'MDSC (Jiang)', xlab = 'Immune classification', out.pch = NA)
  beeswarm(MDSC_jiang ~ immune_cluster, data = patient.summary.region, pch = 16, add = TRUE, col = 'black', method = 'swarm', corral = 'wrap')
  p.val <- wilcox.test(MDSC_jiang ~ immune_cluster, data = patient.summary.region)$p.value
  legend(x = 'bottomright', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
}

# Figure S5D ###
{
  patient.summary.region$tmp_class <- ifelse(patient.summary.region$orig_immune_cluster == 'high' & patient.summary.region$pathology_TILs > 70, yes = 'high_high', no = ifelse(patient.summary.region$orig_immune_cluster == 'high' & patient.summary.region$pathology_TILs <= 20, yes = 'high_low', no = ifelse(patient.summary.region$orig_immune_cluster == 'low' & patient.summary.region$pathology_TILs <= 15, yes = 'low_low', no = ifelse(patient.summary.region$orig_immune_cluster == 'low' & patient.summary.region$pathology_TILs >= 48, yes = 'low_high', no = NA))))
  patient.summary$low_high_region <- ifelse(patient.summary$CRUKid %in% patient.summary.region$CRUKid[which(patient.summary.region$tmp_class == 'low_high')], yes = TRUE, no = FALSE)
  
  par(mfrow = c(1,2))
  
  boxplot(cd8.score.danaher ~ tmp_class, data = patient.summary.region[which(patient.summary.region$tmp_class %in% c('low_low', 'low_high') & patient.summary.region$simple.hist == 'LUAD'),], las = 1, ylab = 'CD8+ score (Danaher)', ylim = c(0,4.5))
  beeswarm(cd8.score.danaher ~ tmp_class, data = patient.summary.region[which(patient.summary.region$tmp_class %in% c('low_low', 'low_high') & patient.summary.region$simple.hist == 'LUAD'),], add=TRUE, pch = 16, cex = 1.4)
  p.val <- wilcox.test(cd8.score.danaher ~ tmp_class, data = patient.summary.region[which(patient.summary.region$tmp_class %in% c('low_low', 'low_high') & patient.summary.region$simple.hist == 'LUAD'),])$p.value
  legend(x = 'topleft', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
  
  boxplot(treg.score.danaher ~ tmp_class, data = patient.summary.region[which(patient.summary.region$tmp_class %in% c('low_low', 'low_high') & patient.summary.region$simple.hist == 'LUAD'),], las = 1, ylab = 'Treg score (Danaher)', ylim = c(0,2.5))
  beeswarm(treg.score.danaher ~ tmp_class, data = patient.summary.region[which(patient.summary.region$tmp_class %in% c('low_low', 'low_high') & patient.summary.region$simple.hist == 'LUAD'),], add=TRUE, pch = 16, cex = 1.4)
  p.val <- wilcox.test(treg.score.danaher ~ tmp_class, data = patient.summary.region[which(patient.summary.region$tmp_class %in% c('low_low', 'low_high') & patient.summary.region$simple.hist == 'LUAD'),])$p.value
  legend(x = 'topleft', legend = paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
}

# Figure S5E
{
  layout(mat = matrix(c(1,2,0,0,3,4), nrow = 2), heights = c(0.9,0.1), widths = c(1,0.2,1))
  par(oma = c(3.5,2,0,1)+0.1, mar = c(0,2.5,0,0)+0.5)
  
  tmp  <- patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),]
  tmp  <- tmp[!is.na(tmp$pathology_TILs),]
  boop <- with(tmp, reorder(CRUKid, pathology_TILs, median, na.rm = TRUE))
  bp   <- boxplot(pathology_TILs ~ boop, ylim = c(0,100), data = tmp, las = 2, space = 0.3, cex.axis = 1, cex.lab = 0.8, names = FALSE, frame.plot = FALSE, at = 1:length(levels(boop)), outpch = NA)
  beeswarm(pathology_TILs ~ boop, data = tmp, pch = 16, add = TRUE, cex = 1.5, pwcol = ifelse(is.na(tmp$orig_immune_cluster), yes = 'light grey', no = ifelse(tmp$orig_immune_cluster == 'high', yes = '#BD2131', no = '#30ABDF')))
  
  mtext(text = 'pathology TILs', side = 2, line = 2.5, cex = 0.8)
  legend('topleft', legend = c('immune high cluster', 'immune low cluster'), pch = 16, col = c('#BD2131','#30ABDF'), bty = 'n', cex = 1)
  legend('bottomright', legend = 'LUAD', bty = 'n')
  beep <- patient.summary$orig_immune_classification[match(levels(boop), patient.summary$CRUKid)]
  bp   <- barplot(rep(1, length(beep)), col = ifelse(is.na(beep), yes = 'darkgrey', no = ifelse(beep == 'low', yes = '#30ABDF', no = ifelse(beep == 'mixed', yes = '#df9b2f', no = '#BD2131'))), space = 0.3, axes = FALSE, border = NA)
  mtext('immune\nclass', side = 2, line = 1, cex = 0.8, las = 1)
  mtext(text = levels(boop), side = 1, line = 1, las = 2, at = bp[,1], cex = 0.7)
  
  
  tmp  <- patient.summary.region[which(patient.summary.region$simple.hist == 'LUSC'),]
  tmp  <- tmp[!is.na(tmp$pathology_TILs),]
  boop <- with(tmp, reorder(CRUKid, pathology_TILs, median, na.rm = TRUE))
  bp   <- boxplot(pathology_TILs ~ boop, data = tmp, ylim = c(0,100), las = 2, space = 0.3, cex.axis = 1, cex.lab = 0.8, names = FALSE, frame.plot = FALSE, at = 1:length(levels(boop)), outpch = NA)
  beeswarm(pathology_TILs ~ boop, data = tmp, pch = 16, add = TRUE, cex = 1.5, pwcol = ifelse(is.na(tmp$immune_cluster), yes = 'light grey', no = ifelse(tmp$immune_cluster == 'high', yes = '#BD2131', no = '#30ABDF')))
  mtext(text = 'pathology TILs', side = 2, line = 2.5, cex = 0.8)
  legend('topleft', legend = c('immune high cluster', 'immune low cluster'), pch = 16, col = c('#BD2131','#30ABDF'), bty = 'n', cex = 1)
  legend('bottomright', legend = 'LUSC', bty = 'n')
  beep <- patient.summary$orig_immune_classification[match(levels(boop), patient.summary$CRUKid)]
  bp   <- barplot(rep(1, length(beep)), col = ifelse(is.na(beep), yes = 'darkgrey', no = ifelse(beep == 'low', yes = '#30ABDF', no = ifelse(beep == 'mixed', yes = '#df9b2f', no = '#BD2131'))), space = 0.3, axes = FALSE, border = NA)
  mtext('immune\nclass', side = 2, line = 1, cex = 0.8, las = 1)
  mtext(text = levels(boop), side = 1, line = 1, las = 2, at = bp[,1], cex = 0.7)

  par(orig_par)
}

# Figure S5F
{
  par(mfrow = c(1,2))
  tmp <- patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),]
  tmp$regionalPattern[which(tmp$regionalPattern %in% c('MP', 'scar only', 'no tumour'))] <- 'other'
  tmp$regionalPattern <- factor(tmp$regionalPattern, levels = c('cribriform', 'solid', 'lepidic', 'other', 'papillary', 'acinar'))
  het_luad_patients <- unique(tmp$CRUKid[which(tmp$orig_immune_classification == 'mixed')])
  het_tmp <- tmp[which(tmp$CRUKid %in% het_luad_patients),c('CRUKid', 'immune_cluster', 'regionalPattern')]
  het_tmp$immune_cluster  <- as.character(het_tmp$immune_cluster)
  het_tmp$regionalPattern <- as.character(het_tmp$regionalPattern)
  het_tmp$regionalPattern[is.na(het_tmp$regionalPattern)] <- 'NA'
  het_tmp <- het_tmp[!is.na(het_tmp$immune_cluster),]
  
  patient.summary$pathologic_heterogeneity <- sapply(patient.summary$CRUKid, FUN = function(x){if(patient.summary$simple.hist[which(patient.summary$CRUKid == x)] %in% c('LUSC', 'other')){return(NA)}; tmp<-patient.summary.region[which(patient.summary.region$CRUKid == x),]; tmp <- tmp[!is.na(tmp$regionalPattern),]; tmp$regionalPattern[which(tmp$regionalPattern %in% c('MP', 'scar only', 'no tumour'))] <- 'other';  if(all(is.na(tmp$regionalPattern))){return(NA)}; ifelse(length(unique(tmp$regionalPattern)) == 1, yes = return('homo'), no = return('het'))})
  p.val <- round(fisher.test(table(patient.summary$pathologic_heterogeneity, patient.summary$orig_immune_classification == 'mixed'))$p.value, digits = 3)
  
  layout(mat = matrix(c(1,2), nrow = 1), widths = c(0.4,0.6))
  
  barplot(table(patient.summary$pathologic_heterogeneity, patient.summary$immune_classification == 'mixed'), names.arg = c('Homogeneous immune', 'Heterogeneous immune'), col = rev(c('#00441b', '#a1d99b')), las = 1, xlab = 'Lung adeno. patients', ylab = 'Number patients', cex.axis = 1, cex.names = 1, cex.lab = 1)
  legend('topright',legend = c(paste('p: ', p.val, sep = ''), 'Heterogeneous pathology', 'Homogeneous pathology'), fill = rev(c('#00441b', '#a1d99b', NA)), bty = 'n', border = NA, cex = 0.8)
  
  plot(x = 0, axes = FALSE, pch = NA, xlim = c(1,30), ylim = c(0,8), ylab = 'Tumor Regions', xlab = '', main = 'Histologies of regions from heterogeneous patients', cex.main = 0.9, frame.plot = FALSE, las =1, cex.lab = 0.9)
  histcols <- setNames(c('#a6cee3','#33a02c', '#ff7f00', '#6a3d9a', '#fb9a99', 'lightgrey', '#a65628'), c('acinar','solid', 'lepidic', 'cribriform', 'papillary', 'NA', 'other'))
  
  for(pat in sort(unique(het_tmp$CRUKid))){
    x = which(sort(unique(het_tmp$CRUKid)) == pat)
    tmp_plot <- het_tmp[which(het_tmp$CRUKid == pat),]
    tmp_plot <- tmp_plot[order(tmp_plot$regionalPattern),]
    tmp_plot <- tmp_plot[order(tmp_plot$immune_cluster, decreasing = TRUE),]
    tmp_plot$immune_color <- ifelse(tmp_plot$immune_cluster == 'low', yes = '#30ABDF', no = '#BD2131')
    tmp_plot$hist_color   <- histcols[match(tmp_plot$regionalPattern, names(histcols))]
    for(i in 1:nrow(tmp_plot)){
      rect(xleft = 2*x-0.8, xright = 2*x, ybottom = i-1, ytop = i, col = tmp_plot$immune_color[i], border = NA)
      rect(xleft = 2*x, xright = 2*x+0.8, ybottom = i-1, ytop = i, col = tmp_plot$hist_color[i], border = NA)
    }
  }
  
  axis(side = 1, at = seq(from = 2, to = 16, by = 2), labels = sort(unique(het_tmp$CRUKid)), las = 2, cex.axis = 0.8)
  axis(side = 2, at = c(0,4,8), labels = c(0,4,8), las = 1, cex.axis = 0.8)
  
  legend(x = 'topleft', legend = c('low', 'high', names(histcols)), fill = c('#30ABDF', '#BD2131', histcols), border = NA, bty = 'n', ncol = 5, cex = 0.8)
}

# Figure S5G
{
  par(mfrow = c(1,2))
  barplot(table(patient.summary$orig_immune_classification)[c('low', 'mixed', 'high')], col = c('#30ABDF', '#df9b2f', '#BD2131'), border = NA, las = 1, ylab = 'Number patients', ylim = c(0,40))
  barplot(table(patient.summary$immune_classification)[c('low', 'mixed', 'high')], col = c('#30ABDF', '#df9b2f', '#BD2131'), border = NA, las = 1, ylab = 'Number patients', ylim = c(0,40))
  
}

# Figure S5H
{
  par(mfrow = c(1,1))
  het_tmb_patients <- patient.summary$CRUKid[which(patient.summary$TMB_het == TRUE)]
  het_tmb_dataframe <- data.frame(row.names = het_tmb_patients)
  het_tmb_dataframe$low_tmb_purity <- sapply(rownames(het_tmb_dataframe), FUN = function(x){
    regions_low_tmb <- patient.summary.region$CRUKidRegion[which(patient.summary.region$CRUKid == x & patient.summary.region$TotalMut_Varscan/30 <= 10)]
    summary(patient.summary.region$purity[which(patient.summary.region$CRUKidRegion  %in% regions_low_tmb)])['Median']
  })
  het_tmb_dataframe$high_tmb_purity <- sapply(rownames(het_tmb_dataframe), FUN = function(x){
    regions_high_tmb <- patient.summary.region$CRUKidRegion[which(patient.summary.region$CRUKid == x & patient.summary.region$TotalMut_Varscan/30 > 10)]
    summary(patient.summary.region$purity[which(patient.summary.region$CRUKidRegion  %in% regions_high_tmb)])['Median']
  })
  
  p.val <- t.test(het_tmb_dataframe$low_tmb_purity, het_tmb_dataframe$high_tmb_purity, paired = TRUE)$p.value
  
  stripchart(het_tmb_dataframe, vertical = TRUE, pch= 16, at = c(0.9, 1.2), xlim =c(0.7,1.4), las = 1, group.names = c('Low TMB regions', 'High TMB regions'), cex = 1.5, ylab = 'Tumor region purity')
  segments(x0 = 0.9, y0 = het_tmb_dataframe$low_tmb_purity, x1 = 1.2, y1 = het_tmb_dataframe$high_tmb_purity, col = '#df9b2f', lwd = 3)
  legend(legend = paste('paired t-test p: ', round(p.val, digits = 2), sep = ''), x = 'topleft', bty = 'n', cex = 0.9)
}

# Figure S6A
{
  layout(mat = matrix(c(1,2,3,0), nrow = 4), heights = c(0.9,0.1,0.1, 0.1), widths = c(1))
  par(oma = c(3.5,2,0,1)+0.1, mar = c(0,2.5,0,0)+0.5)
  
  tmp  <- patient.summary.region[which(patient.summary.region$simple.hist %in% c('LUAD', 'LUSC')),]
  tmp  <- tmp[!is.na(tmp$TIDE_score),]
  tmp  <- tmp[which(tmp$number.regions.rna.accepted > 1),]
  boop <- with(tmp, reorder(CRUKid, TIDE_score, FUN = function(x) -min(x, na.rm = TRUE)))
  bp   <- boxplot(TIDE_score ~ boop, data = tmp, ylim = c(-2,3), las = 2, space = 0.3, cex.axis = 1, cex.lab = 0.8, names = FALSE, frame.plot = FALSE, at = 1:length(levels(boop)), outpch = NA)
  beeswarm(TIDE_score ~ boop, data = tmp, pch = 16, add = TRUE, cex = 1, pwcol = ifelse(is.na(tmp$immune_cluster), yes = 'light grey', no = ifelse(tmp$immune_cluster == 'high', yes = '#BD2131', no = '#30ABDF')))
  abline(h = 0, col = 'black', lty = 2)
  text(x = 0.1, y = 0.1, labels = 'TIDE threshold', pos = 4, cex = 0.7)
  mtext(text = 'TIDE Score', side = 2, line = 2.5, cex = 1)
  legend('topleft', legend = c('immune high cluster', 'immune low cluster'), pch = 16, col = c('#BD2131','#30ABDF'), bty = 'n', cex = 1)
  #legend('bottomright', legend = cancer, bty = 'n')
  
  beep <- patient.summary$immune_classification[match(levels(boop), patient.summary$CRUKid)]
  bp   <- barplot(rep(1, length(beep)), col = ifelse(beep == 'low', yes = '#30ABDF', no = ifelse(beep == 'mixed', yes = '#df9b2f', no = '#BD2131')), space = 0.3, axes = FALSE, border = NA)
  mtext('immune\nclass', side = 2, line = 1, cex = 0.8, las = 1)
  
  beep       <- patient.summary$TIDE_enrichment_group[match(levels(boop), patient.summary$CRUKid)]
  beep_color <- ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'False', yes = '#e0f3f8', no = ifelse(beep == 'True', yes = '#4575b4', no = '#df9b2f')))
  bp         <- barplot(rep(1, length(beep)), col = beep_color, space = 0.3, axes = FALSE, border = NA)
  mtext('TIDE', side = 2, line = 1, cex = 0.8, las = 1)

  par(orig_par)
}

# Figure S6B
{
  layout(mat = matrix(c(1,2,3,0), nrow = 4), heights = c(0.9,0.1,0.1, 0.1, 0.1), widths = c(1))
  par(oma = c(3.5,2,0,1)+0.1, mar = c(0,2.5,0,0)+0.5)
  
  tmp  <- patient.summary.region[which(patient.summary.region$simple.hist %in% c('LUAD', 'LUSC')),]
  tmp  <- tmp[!is.na(tmp$ipres_enrichment),]
  tmp  <- tmp[which(tmp$number.regions.rna.accepted > 1),]
  boop <- with(tmp, reorder(CRUKid, ipres_enrichment, max, na.rm = TRUE))
  bp   <- boxplot(ipres_enrichment ~ boop, data = tmp, ylim = c(-2,2), las = 2, space = 0.3, cex.axis = 1, cex.lab = 0.8, names = FALSE, frame.plot = FALSE, at = 1:length(levels(boop)), outpch = NA)
  beeswarm(ipres_enrichment ~ boop, data = tmp, pch = 16, add = TRUE, cex = 1.2, pwcol = ifelse(is.na(tmp$immune_cluster), yes = 'light grey', no = ifelse(tmp$immune_cluster == 'high', yes = '#BD2131', no = '#30ABDF')))
  abline(h = 0.35, col = 'black', lty = 2)
  text(x = 0.1, y = 0.4, labels = 'IPRES threshold', pos = 4, cex = 0.7)
  mtext(text = 'IPRES Score', side = 2, line = 2.5, cex = 1)
  legend('topleft', legend = c('immune high cluster', 'immune low cluster'), pch = 16, col = c('#BD2131','#30ABDF'), bty = 'n', cex = 1)
  
  beep       <- patient.summary$ipres_enrichment_group[match(levels(boop), patient.summary$CRUKid)]
  beep_color <- ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'FALSE', yes = '#e0f3f8', no = ifelse(beep == 'TRUE', yes = '#4575b4', no = '#df9b2f')))
  bp         <- barplot(rep(1, length(beep)), col = beep_color, space = 0.3, axes = FALSE, border = NA)
  mtext('IPRES', side = 2, line = 1, cex = 0.8, las = 1)
  
  beep <- patient.summary$immune_classification[match(levels(boop), patient.summary$CRUKid)]
  bp   <- barplot(rep(1, length(beep)), col = ifelse(beep == 'low', yes = '#30ABDF', no = ifelse(beep == 'mixed', yes = '#df9b2f', no = '#BD2131')), space = 0.3, axes = FALSE, border = NA)
  mtext('immune\nclass', side = 2, line = 1, cex = 0.8, las = 1) 
  
  par(orig_par)
}

# Figure S6C-D
{
  layout(mat = matrix(c(1,2,3,3), nrow = 2), heights = c(0.9, 0.1), widths = c(0.7, 0.3))
  par(oma = c(3.5,2,0,1)+0.1, mar = c(0,2.5,0,0)+0.5)
  
  tmp  <- patient.summary.region[which(patient.summary.region$simple.hist %in% c('LUAD', 'LUSC')),]
  tmp  <- tmp[!is.na(tmp$ifng_signature_expanded),]
  tmp  <- tmp[which(tmp$number.regions.rna.accepted > 1),]
  boop <- with(tmp, reorder(CRUKid, ifng_signature_expanded, max, na.rm = TRUE))
  bp   <- boxplot(ifng_signature_expanded ~ boop, data = tmp, ylim = c(4,10), las = 2, space = 0.3, cex.axis = 1, cex.lab = 0.8, names = FALSE, frame.plot = FALSE, at = 1:length(levels(boop)), outpch = NA)
  beeswarm(ifng_signature_expanded ~ boop, data = tmp, pch = 16, add = TRUE, cex = 1.2, pwcol = ifelse(is.na(tmp$immune_cluster), yes = 'light grey', no = ifelse(tmp$immune_cluster == 'high', yes = '#BD2131', no = '#30ABDF')))
  mtext(text = 'IFNg signature, expanded', side = 2, line = 2.5, cex = 1)
  legend('topleft', legend = c('immune high cluster', 'immune low cluster'), pch = 16, col = c('#BD2131','#30ABDF'), bty = 'n', cex = 1)

  beep <- patient.summary$immune_classification[match(levels(boop), patient.summary$CRUKid)]
  bp   <- barplot(rep(1, length(beep)), col = ifelse(beep == 'low', yes = '#30ABDF', no = ifelse(beep == 'mixed', yes = '#df9b2f', no = '#BD2131')), space = 0.3, axes = FALSE, border = NA)
  mtext('immune\nclass', side = 2, line = 1, cex = 0.8, las = 1)
  
  boxplot(Ayers_diff ~ immune_classification == 'mixed', patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),], las = 1, xlab = 'Immune heterogeneous', ylab = 'IFNg signature, expanded (Ayers)')
  beeswarm(Ayers_diff ~ immune_classification == 'mixed', patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),], pch = 16, add=TRUE)
  p.val <- wilcox.test(Ayers_diff ~ immune_classification == 'mixed', patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),])
  legend('topleft', paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), sep = ''), bty = 'n', cex = 1)
  
  par(orig_par)
}

# Figure S6E
{
  layout(mat = matrix(c(1,2), nrow = 2), heights = c(0.9,0.1), widths = c(1))
  par(oma = c(3.5,2,0,1)+0.1, mar = c(0,2.5,0,0)+0.5)
  
  tmp  <- patient.summary.region[which(patient.summary.region$simple.hist != 'other'),]
  
  boop <- with(tmp, reorder(CRUKid, TotalMut_Varscan, max, na.rm = TRUE))
  bp   <- boxplot(TotalMut_Varscan/30 ~ boop, data = tmp, ylim = c(0.5,100), las = 2, space = 0.3, cex.axis = 1, cex.lab = 0.8, names = FALSE, frame.plot = FALSE, at = 1:length(levels(boop)), outpch = NA, log = 'y', col = 'white')
  stripchart(TotalMut_Varscan/30 ~ boop, data = tmp, pch = 16, add = TRUE, cex = 1, vertical=TRUE)
  abline(h = 10, col = 'black', lty = 2)
  text(x = 0.1, y = 11, labels = 'TMB threshold', pos = 4, cex = 0.7)
  mtext(text = 'TMB (mut/Mb) [VarScan]', side = 2, line = 2.5, cex = 1)
  
  beep <- patient.summary$TMB_het_col[match(levels(boop), patient.summary$CRUKid)]
  bp   <- barplot(rep(1, length(beep)), col = beep, space = 0.3, axes = FALSE, border = NA)
  mtext('TMB +', side = 2, line = 1, cex = 0.8, las = 1)
  
  par(orig_par)
}

# Figure S6F
{
  par(mfrow = c(5,1), mar = c(1,4,1,2))
  
  tmp <- patient.summary[,c('simple.hist', 'immune_classification', 'ipres_enrichment_group', 'TIDE_enrichment_group', 'TMB_any', 'TMB_het')]
  tmp$TMB_group <- ifelse(tmp$TMB_any & tmp$TMB_het, yes = 'mixed', no = ifelse(tmp$TMB_any, yes = 'TRUE', no ='FALSE'))
  tmp$TIDE_enrichment_group <- gsub(pattern = 'True', replacement = 'TRUE', x = tmp$TIDE_enrichment_group)
  tmp$TIDE_enrichment_group <- gsub(pattern = 'False', replacement = 'FALSE', x = tmp$TIDE_enrichment_group)
  tmp <- tmp[which(tmp$simple.hist != 'other'),]
  
  tmp$num_pos <- rowSums(tmp[,c('TMB_group', 'TIDE_enrichment_group', 'ipres_enrichment_group')]=='TRUE', na.rm = TRUE)
  
  tmp <- tmp[with(tmp, order(immune_classification, num_pos, TMB_group, TIDE_enrichment_group, ipres_enrichment_group, decreasing = TRUE)),]
  
  beep       <- tmp$simple.hist
  beep_color <- ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'LUAD', yes = 'dark blue', no = ifelse(beep == 'LUSC', yes = 'dark red', no = 'lightgrey')))
  bp         <- barplot(rep(1, length(beep)), col = beep_color, space = 0.3, axes = FALSE, border = NA)
  mtext('Hist.', side = 2, line = 1, cex = 0.8, las = 1)
  
  beep <- tmp$immune_classification
  bp   <- barplot(rep(1, length(beep)), col = ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'low', yes = '#30ABDF', no = ifelse(beep == 'mixed', yes = '#df9b2f', no = '#BD2131'))), space = 0.3, axes = FALSE, border = NA)
  mtext('immune\nclass', side = 2, line = 1, cex = 0.8, las = 1)
  
  beep       <- tmp$TMB_group
  beep_color <- ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'FALSE', yes = '#e0f3f8', no = ifelse(beep == 'TRUE', yes = '#4575b4', no = '#df9b2f')))
  bp         <- barplot(rep(1, length(beep)), col = beep_color, space = 0.3, axes = FALSE, border = NA)
  mtext('TMB+', side = 2, line = 1, cex = 0.8, las = 1)
  
  beep       <- tmp$TIDE_enrichment_group
  beep_color <- ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'FALSE', yes = '#e0f3f8', no = ifelse(beep == 'TRUE', yes = '#4575b4', no = '#df9b2f')))
  bp         <- barplot(rep(1, length(beep)), col = beep_color, space = 0.3, axes = FALSE, border = NA)
  mtext('TIDE', side = 2, line = 1, cex = 0.8, las = 1)
  
  beep       <- tmp$ipres_enrichment_group
  beep_color <- ifelse(is.na(beep), yes = 'lightgrey', no = ifelse(beep == 'FALSE', yes = '#e0f3f8', no = ifelse(beep == 'TRUE', yes = '#4575b4', no = '#df9b2f')))
  bp         <- barplot(rep(1, length(beep)), col = beep_color, space = 0.3, axes = FALSE, border = NA)
  mtext('IPRES', side = 2, line = 1, cex = 0.8, las = 1) 
  
  par(orig_par)
}

# Figure S7A ###
{
  par(mfrow = c(1,1))
  plot(total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'cn'], total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'immune'], xlab = 'Pairwise cn distance', ylab = 'Pairwise immune distance', pch = 16, col = total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'color'], las = 1)
  p.val <- cor.test(total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'cn'], total_comparisons[which(total_comparisons$simple.hist == 'LUAD'),'immune'], method = 's')
  legend(x = 'topleft', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\n', 'rho: ', formatC(p.val$estimate, digits = 1, format = 'e'), sep= ''), bty = 'n')
  
  points(total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'cn'], total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'immune'], xlab = 'Pairwise cn distance', ylab = 'Pairwise immune distance', pch = 16, col = total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'color'])
  p.val <- cor.test(total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'cn'], total_comparisons[which(total_comparisons$simple.hist == 'LUSC'),'immune'], method = 's')
  legend(x = 'topright', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\n', 'rho: ', formatC(p.val$estimate, digits = 1, format = 'e'), sep= ''), bty = 'n')  
  
}

# Figure S7C
{
  par(mfrow = c(1,2))
  boxplot((nk.score.danaher)/total.til.score.danahaer ~ !homo_c, data = patient.summary.region[which(patient.summary.region$hla_c_lost == TRUE),], las = 1, ylab = 'NK cell', xlab = 'C1/C2 heterozygosity', main = 'Patients with HLA-C LOH', ylim = c(0,1))
  beeswarm((nk.score.danaher)/total.til.score.danahaer ~ !homo_c, data = patient.summary.region[which(patient.summary.region$hla_c_lost == TRUE),], add=TRUE, pch = 16, col = 'darkgreen')
  p.val <- wilcox.test((nk.score.danaher/total.til.score.danahaer) ~ homo_c, data = patient.summary.region[which(patient.summary.region$hla_c_lost == TRUE),])$p.value
  legend('topleft', paste('p: ', formatC(p.val, digits = 1, format = 'e'), sep = ''), bty = 'n', cex = 0.8)
  legend('topright', 'tracerx', cex = 0.8, bty = 'n')
  
  boxplot((nk.score.danaher)/total.til.score.danahaer ~ !homo_c, data = patient.summary.region[which(patient.summary.region$hla_c_lost == FALSE),], las = 1, ylab = 'NK cell', xlab = 'C1/C2 heterozygosity', main = 'Patients without HLA-C LOH', ylim = c(0,1))
  beeswarm((nk.score.danaher)/total.til.score.danahaer ~ !homo_c, data = patient.summary.region[which(patient.summary.region$hla_c_lost == FALSE),], add=TRUE, pch = 16, col = 'darkgreen')
  p.val <- wilcox.test((nk.score.danaher)/total.til.score.danahaer ~ homo_c, data = patient.summary.region[which(patient.summary.region$hla_c_lost == FALSE),])$p.value
  legend('topleft', paste('p: ', formatC(p.val, digits = 1, format = 'e'), sep = ''), bty = 'n', cex = 0.8)
  
}

# Figure S7D-E
{
  par(mfrow = c(1,2))
  plot.with.confidence(patient.summary.region$shannon.ith[which(patient.summary.region$simple.hist == 'LUAD')], patient.summary.region$cd8.score.danaher[which(patient.summary.region$simple.hist == 'LUAD')], pch = 16, xlab = 'Shannon ITH', ylab = 'CD8+ Danaher', main = 'Lung adeno.', col = 'dark blue', cex = 1.6)
  plot.with.confidence(patient.summary.region$shannon.ith[which(patient.summary.region$simple.hist == 'LUSC')], patient.summary.region$cd8.score.danaher[which(patient.summary.region$simple.hist == 'LUSC')], pch = 16, xlab = 'Shannon ITH', ylab = 'CD8+ Danaher', main = 'Lung squam.', col = 'dark red', cex = 1.6)
}

# Figure S7F
{
  par(mfrow = c(1,1))
  plot(patient.summary.region$purity[which(patient.summary.region$simple.hist == 'LUAD')], patient.summary.region$pathology_TILs[which(patient.summary.region$simple.hist == 'LUAD')], pch = 16, col = 'dark blue', ylab = 'pathology estimated TILs', xlab = 'tumor purity', las = 1, main = '', cex.main = 1, font.main = 1, xlim = c(0,1))
  p.val <- cor.test(patient.summary.region$purity[which(patient.summary.region$simple.hist == 'LUAD')], patient.summary.region$pathology_TILs[which(patient.summary.region$simple.hist == 'LUAD')], method = 's')
  legend(x = 'topright', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\n', 'rho: ', formatC(p.val$estimate, digits = 1, format = 'e'), sep= ''), bty = 'n', cex = 0.7)

  points(patient.summary.region$purity[which(patient.summary.region$simple.hist == 'LUSC')], patient.summary.region$pathology_TILs[which(patient.summary.region$simple.hist == 'LUSC')], pch = 16, col = 'dark red', ylab = 'pathology estimated TILs', xlab = 'tumor purity')
  p.val <- cor.test(patient.summary.region$purity[which(patient.summary.region$simple.hist == 'LUSC')], patient.summary.region$pathology_TILs[which(patient.summary.region$simple.hist == 'LUSC')], method = 's')
  legend(x = 'bottomright', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\n', 'rho: ', formatC(p.val$estimate, digits = 1, format = 'e'), sep= ''), bty = 'n', cex = 0.7)
}

# Figure S7G
{
  par(mfrow = c(1,1))
  tmp <- gsub(pattern = '0%', replacement = '#30ABDF', x = gsub(pattern = '50%', replacement = '#BD2131', x = patient.summary.region$TIL_quart[which(patient.summary.region$simple.hist == 'LUAD')]))
  boxplot(shannon.ith ~ TIL_cat, data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),], las = 1, main = 'Lung adeno.')
  beeswarm(shannon.ith ~ TIL_cat, data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),], add = TRUE, pch = 16, pwcol = tmp, cex = 1.6)
  p.val <- wilcox.test(shannon.ith ~ TIL_cat == 'low', data = patient.summary.region[which(patient.summary.region$simple.hist == 'LUAD'),])$p.value
  legend('bottomright', paste('p -- low-mixed/high : ', round(p.val, digits = 2), sep = ''), bty = 'n')
}

# Figure S7H-I
{
  tmp <- patient.summary[,c('ie.all.updated', 'ie.clonal.updated', 'ie.subclonal.updated', 'simple.hist', 'immune_classification', 'anyCPN_remove')]
  tmp <- tmp[complete.cases(tmp),]
  
  par(mfrow = c(1,2))
  boxplot(ie.all.updated ~ as.character(simple.hist), data = tmp[which(tmp$simple.hist %in% c('LUAD', 'LUSC')),], las = 1, ylab = 'Immunoediting score')
  beeswarm(ie.all.updated ~ as.character(simple.hist), data = tmp[which(tmp$simple.hist %in% c('LUAD', 'LUSC')),], pwcol = ifelse(tmp$simple.hist[which(tmp$simple.hist %in% c('LUAD', 'LUSC'))] == 'LUAD', yes = 'dark blue', no = 'dark red'), pch = 16, add =TRUE)
  p.val <- wilcox.test(ie.all.updated ~ as.character(simple.hist), data = tmp[which(tmp$simple.hist %in% c('LUAD', 'LUSC')),])$p.value
  legend('bottomright', paste('p: ', formatC(p.val, format = 'e', digits = 1), sep = ''), bty = 'n', cex = 1.3)
  
  p.val <- cor.test(as.numeric(patient.summary$ie.all.updated[which(patient.summary$simple.hist != 'other')]),  as.numeric(patient.summary$number.unique.hlas[which(patient.summary$simple.hist != 'other')]), method= 's')
  boxplot(ie.all.updated ~ number.unique.hlas, data = patient.summary[which(patient.summary$simple.hist != 'other'),], las = 1, ylab = 'Immunoediting score', xlab = 'Number unique HLAs')
  beeswarm(ie.all.updated ~ number.unique.hlas, data = patient.summary[which(patient.summary$simple.hist != 'other'),], pch = 16, add =TRUE)
  legend(x = 'topleft', legend = paste('p: ', formatC(p.val$p.value, digits = 1, format = 'e'), '\nrho: ', round(p.val$estimate, digits = 2), sep = ''), bty = 'n', cex = 1.3)
}

# Figure S7J
{
  
  tmp <- mutTableAll.together.filtered[which(mutTableAll.together.filtered$ITHState_RNA %in% c(1,2,3, 'FALSE')),] 
  
  {
    par(mfrow = c(1,1))
    ests   <- c()
    p.vals <- c()
    ors1    <- c()
    ors2    <- c()
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == TRUE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$any.HLA.loss == FALSE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == TRUE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('high') & tmp$any.HLA.loss == FALSE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == TRUE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('mixed') & tmp$any.HLA.loss == FALSE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == TRUE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == TRUE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    tmp2 <- table(tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == FALSE),'lost_in_rna'], tmp[which(tmp$nonsynonymous == TRUE & tmp$immune_classification %in% c('low') & tmp$any.HLA.loss == FALSE),'strong.count'])
    f.test  <- fisher.test(tmp2)
    ests    <- c(ests, f.test$estimate)
    p.vals  <- c(p.vals, f.test$p.value)
    ors1    <- c(ors1, f.test$conf.int[1])
    ors2    <- c(ors2, f.test$conf.int[2])
    
    ests <- 1/ests
    ors1 <- 1/ors1
    ors2 <- 1/ors2
    
    bp <- barplot(ests, col = c(sapply(X = c('#00441b', '#BD2131', '#df9b2f', '#30ABDF'), FUN = function(x) rep(x, 2))), border = FALSE, las = 1, cex.axis = 0.9, cex.names = 0.7, ylab = 'OR -- neoantigen expression lost in RNA', ylim = c(0,1.5), names = '', axes = FALSE, space = rep(c(0.6,0.2),4))
    axis(side = 2, at = c(0,0.5, 1.0, 1.5), labels = c(0, 0.5, 1.0, 1.5), las = 1, cex.axis = 0.8)
    segments(x0 = bp[1,1], x1 = bp[1,1], y0 = ors1[1], y1 = ors2[1])
    segments(x0 = bp[2,1], x1 = bp[2,1], y0 = ors1[2], y1 = ors2[2])
    segments(x0 = bp[3,1], x1 = bp[3,1], y0 = ors1[3], y1 = ors2[3])
    segments(x0 = bp[4,1], x1 = bp[4,1], y0 = ors1[4], y1 = ors2[4])
    segments(x0 = bp[5,1], x1 = bp[5,1], y0 = ors1[5], y1 = ors2[5])
    segments(x0 = bp[6,1], x1 = bp[6,1], y0 = ors1[6], y1 = ors2[6])
    segments(x0 = bp[7,1], x1 = bp[7,1], y0 = ors1[7], y1 = ors2[7])
    segments(x0 = bp[8,1], x1 = bp[8,1], y0 = ors1[8], y1 = ors2[8])
    
    mtext(text = paste('p = ', round(p.vals, digits = 3), sep = ''), side = 1, line = 2, cex = 0.7, at = bp[,1])
    mtext(text = c('+HLA LOH', '-HLA LOH'), side = 1, line = 1, cex = 0.7, at = bp[,1])
    tmp_names <- factor(rep(c('All', 'High', 'Hetero.', 'Low'), each = 2), levels = c('All', 'High', 'Hetero.', 'Low'))
    mtext(text = c('All', 'High', 'Hetero.', 'Low'), side = 1, line = 3, cex = 0.7, at = aggregate(bp, by = list(tmp_names), mean)$V1)
    abline(h = 1, lty = 2, col = 'black')
    
  }
}

# Figure S7K
{
  par(mfrow = c(1,2))
  neo <- table(mutTableRegion$mutant.expressed.binary[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$count == 1)], mutTableRegion$expression.tpm[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$count == 1)] > 1)
  non_neo <- table(mutTableRegion$mutant.expressed.binary[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$count == 0)], mutTableRegion$expression.tpm[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$count == 0)] > 1)
  
  k <- rbind(rev(non_neo[1,]), rev(neo[1,]))
  rownames(k) <- c('non-neo', 'neo')
  
  f.test <- fisher.test(rbind(rev(non_neo[1,]), rev(neo[1,])))
  barplot(t(t(k)/rowSums(t(k))), names.arg = c('Gene expressed', 'Gene not expressed'), las = 1, col = c('#a6bddb', '#2b8cbe'), main = 'Neo')
  legend(x = 'bottomleft', fill =c('#a6bddb', '#2b8cbe'), border = NA, bty = 'n', legend = c('Non-neo', 'Neo'), xpd = TRUE, inset = c(0,-0.35))
  legend(x = 'bottomleft', border = NA, bty = 'n', legend = paste('p: ', formatC(f.test$p.value, digits = 1, format = 'e'), '\nOR: ', round(f.test$estimate, digits = 1), sep = ''), xpd = FALSE)
  
  # strong neo
  neo <- table(mutTableRegion$mutant.expressed.binary[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$strong.count == 1)], mutTableRegion$expression.tpm[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$strong.count == 1)] > 1)
  non_neo <- table(mutTableRegion$mutant.expressed.binary[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$strong.count == 0)], mutTableRegion$expression.tpm[which(mutTableRegion$nonsynonymous == TRUE & mutTableRegion$strong.count == 0)] > 1)
  
  k <- rbind(rev(non_neo[1,]), rev(neo[1,]))
  rownames(k) <- c('non-neo', 'neo')
  
  f.test <- fisher.test(rbind(rev(non_neo[1,]), rev(neo[1,])))
  barplot(t(t(k)/rowSums(t(k))), names.arg = c('Gene expressed', 'Gene not expressed'), las = 1, col = c('#a6bddb', '#2b8cbe'), main = 'Strong neo')
  legend(x = 'bottomleft', fill =c('#a6bddb', '#2b8cbe'), border = NA, bty = 'n', legend = c('Non-neo', 'Neo'), xpd = TRUE, inset = c(0,-0.35))
  legend(x = 'bottomleft', border = NA, bty = 'n', legend = paste('p: ', formatC(f.test$p.value, digits = 1, format = 'e'), '\nOR: ', round(f.test$estimate, digits = 1), sep = ''), xpd = FALSE)
}

# Figure S8
{
  for(cancer in c('LUAD', 'LUSC')){
    if(cancer == 'LUAD'){
      colOrder <- c("CRUK0020", "CRUK0016", "CRUK0039", "CRUK0051", "CRUK0027", "CRUK0060", "CRUK0024", "CRUK0035", "CRUK0047", "CRUK0052", "CRUK0026", "CRUK0008", "CRUK0038", "CRUK0045", "CRUK0031", "CRUK0034", "CRUK0006", "CRUK0033", "CRUK0030", "CRUK0015", "CRUK0014", "CRUK0019", "CRUK0007", "CRUK0040", "CRUK0017", "CRUK0001", "CRUK0009", "CRUK0003", "CRUK0029", "CRUK0061", "CRUK0044", "CRUK0004", "CRUK0012", "CRUK0005", "CRUK0041", "CRUK0018", "CRUK0046", "CRUK0048", "CRUK0032", "CRUK0010", "CRUK0013", "CRUK0028", "CRUK0002", "CRUK0022", "CRUK0037", "CRUK0055", "CRUK0023", "CRUK0036", "CRUK0042", "CRUK0053", "CRUK0057", "CRUK0025", "CRUK0050", "CRUK0021", "CRUK0058", "CRUK0011", "CRUK0043", "CRUK0049", "CRUK0054", "CRUK0056", "CRUK0059")
    }
    if(cancer == 'LUSC'){
      colOrder <- c(c("CRUK0086","CRUK0074","CRUK0068","CRUK0079","CRUK0069","CRUK0064","CRUK0083","CRUK0072","CRUK0073","CRUK0081","CRUK0065","CRUK0078","CRUK0067","CRUK0070","CRUK0062","CRUK0082","CRUK0076","CRUK0071","CRUK0066","CRUK0084","CRUK0090","CRUK0075","CRUK0063","CRUK0089","CRUK0085","CRUK0088","CRUK0087","CRUK0092","CRUK0077","CRUK0080","CRUK0093","CRUK0091"))
    }
    
    tmp_df <- data.frame(colnames(all_uncondensed_alterations)[which(substr(colnames(all_uncondensed_alterations), 1, 8) %in% colOrder)], stringsAsFactors = FALSE)
    tmp_df$SampleID <- substr(tmp_df[,1], 1, 8)
    tmp_df$SampleID <- factor(tmp_df$SampleID, levels = colOrder)
    tmp_df <- tmp_df[order(tmp_df$SampleID),]
    
    uncondensed_alterations <- all_uncondensed_alterations[,tmp_df[,1]]
    
    # mutations annotated as 0.5
    uncondensed_beep <- uncondensed_alterations
    uncondensed_beep[which(uncondensed_beep == -1)] <- 0
    uncondensed_beep[which(uncondensed_beep == 1)] <- 0
    
    # plot
    layout(matrix(1:(nrow(uncondensed_alterations)), ncol = 1))
    par(mar = c(2,2,1,2), oma = c(4,4,4,4))
    tmp_acceptedRegions_df <- acceptedRegions_df[match(colnames(uncondensed_alterations), acceptedRegions_df$CRUKidRegion),]
    
    sapply(1:nrow(uncondensed_alterations), FUN = function(x){
      bp <-  barplot(rep(1,ncol(uncondensed_alterations)),border = NA,axes = FALSE,
                     col = ifelse(as.numeric(uncondensed_alterations[x,]) == 1, yes = 'red', no = ifelse(as.numeric(uncondensed_alterations[x,]) < 0, yes = 'darkblue', no = 'lightgrey')), 
                     space = unlist(sapply(tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)], FUN = function(x){c(2,rep(0.2,(x-1)))}))
      )
      rect(xleft = bp[,1]-0.15, xright = bp[,1]+0.15, ybottom = 0.2, ytop = uncondensed_beep[x,]*1 + 0.2,
           col = '#008000', border = NA)
      mtext(text = rownames(uncondensed_alterations)[x], side=2,line = -25, cex = 5, las = 2, outer = FALSE, col = 'black')
    })
    bp <-  barplot(rep(1,ncol(uncondensed_alterations)),border = NA,axes = FALSE,space = unlist(sapply(tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)], FUN = function(x){c(2,rep(0.2,(x-1)))})), plot = FALSE)
    nr          <- tmp_acceptedRegions_df$nregions[!duplicated(tmp_acceptedRegions_df$CRUKid)]
    names(nr)   <- unique(tmp_acceptedRegions_df$CRUKid)
  }
  par(orig_par)
}



# Figure S9A,C,E
{
  par(mfrow = c(3,2))
  tmp     <- patient.summary[which(patient.summary$simple.hist == 'LUAD'),]
  tmp$tmp <- ifelse(tmp$ClonalNeo > summary(tmp$ClonalNeo)['3rd Qu.'], yes = 'high', no = 'low')
  
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ tmp, 
           data = tmp, 
           main = "LUAD -- survival by upper quartile clonal neo", 
           xlab = 'Time (days)', 
           ylab = "PFS",
           col = c('darkblue', 'darkred'),
           las = 1,
           legend.pos = NA,
           lwd = 1)
  
  tmp$tmp <- ifelse(tmp$SubclonalNeo > summary(tmp$SubclonalNeo)['3rd Qu.'], yes = 'high', no = 'low')
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ tmp, 
           data = tmp, 
           main = "LUAD -- survival by upper quartile subclonal neo", 
           xlab = 'Time (days)', 
           ylab = "PFS", 
           col = c('darkblue', 'darkred'),
           las = 1,
           legend.pos = NA,
           lwd = 1)
  
  
  tmp$tmp <- ifelse(tmp$TotalNeo > summary(tmp$TotalNeo)['3rd Qu.'], yes = 'high', no = 'low')
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ tmp, 
           data = tmp, 
           main = "LUAD -- survival by upper quartile total neo", 
           xlab = 'Time (days)', 
           ylab = "PFS",
           col = c('darkblue', 'darkred'),
           las = 1,
           legend.pos = NA,
           lwd = 1)
  
  tmp     <- patient.summary[which(patient.summary$simple.hist == 'LUSC'),]
  tmp$tmp <- ifelse(tmp$ClonalNeo > summary(tmp$ClonalNeo)['3rd Qu.'], yes = 'high', no = 'low')
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ tmp, 
           data = tmp, 
           main = "LUSC -- survival by upper quartile clonal neo", 
           xlab = 'Time (days)', 
           ylab = "PFS", 
           col = c('darkblue', 'darkred'),
           las = 1,
           legend.pos = NA,
           lwd = 1)
  
  tmp$tmp <- ifelse(tmp$SubclonalNeo > summary(tmp$SubclonalNeo)['3rd Qu.'], yes = 'high', no = 'low')
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ tmp, 
           data = tmp, 
           main = "LUSC -- survival by upper quartile subclonal neo", 
           xlab = 'Time (days)', 
           ylab = "PFS", 
           col = c('darkblue', 'darkred'),
           las = 1,
           legend.pos = NA,
           lwd = 1)
  
  tmp$tmp <- ifelse(tmp$TotalNeo > summary(tmp$TotalNeo)['3rd Qu.'], yes = 'high', no = 'low')
  survplot(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ tmp, 
           data = tmp, 
           main = "LUAD -- survival by upper quartile total neo", 
           xlab = 'Time (days)', 
           ylab = "PFS",
           col = c('darkblue', 'darkred'),
           las = 1,
           legend.pos = NA,
           lwd = 1)
}


# Figure S9B,D ###
{
  par(mfrow = c(2,2), xpd = TRUE)
  
  # # luad
  tmp        <- patient.summary[which(patient.summary$simple.hist == 'LUAD'),]
  luad_pvals <- data.frame(row.names = 1:max(tmp$ClonalNeo))
  luad_pvals$hr <- NA
  luad_pvals$p  <- NA
  luad_pvals$cil <- NA
  luad_pvals$cih <- NA
  luad_pvals$nu <- NA
  luad_pvals$hr_s <- NA
  luad_pvals$p_s  <- NA
  luad_pvals$cil_s <- NA
  luad_pvals$cih_s <- NA
  luad_pvals$nu_s <- NA
  for(n in 1:max(tmp$ClonalNeo)){
    test <- summary(coxph(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ ClonalNeo > n, data = tmp))
    hr  <- test$conf.int[1]
    cil <- test$conf.int[3]
    cih <- test$conf.int[4]
    p   <- test$sctest['pvalue']
    nu  <- survfit(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ ClonalNeo > n, data = tmp)$n[2]
    test <- summary(coxph(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ SubclonalNeo > n, data = tmp))
    hr_s <- test$conf.int[1]
    cil_s <- test$conf.int[3]
    cih_s <- test$conf.int[4]
    p_s  <- test$sctest['pvalue']
    nu_s <- survfit(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ SubclonalNeo > n, data = tmp)$n[2]
    luad_pvals[as.character(n),'hr'] <- hr
    luad_pvals[as.character(n),'p']  <- p
    luad_pvals[as.character(n),'cil']  <- cil
    luad_pvals[as.character(n),'cih']  <- cih
    luad_pvals[as.character(n),'nu'] <- nu
    luad_pvals[as.character(n),'hr_s'] <- hr_s
    luad_pvals[as.character(n),'p_s']  <- p_s
    luad_pvals[as.character(n),'cil_s']  <- cil_s
    luad_pvals[as.character(n),'cih_s']  <- cih_s
    luad_pvals[as.character(n),'nu_s'] <- nu_s
  }
  keep_luad <- luad_pvals
  luad_pvals <- luad_pvals[20:500,]
  
  rbPal           <- colorRampPalette(c('white','darkblue'))
  luad_pvals$ncol <- rbPal(30)[as.numeric(cut(luad_pvals$nu,breaks = 30))]
  luad_pvals$col <- ifelse(luad_pvals$p < 0.05, yes = 'red', no = 'black')
  luad_pvals$cih[which(luad_pvals$cih > 4)] <- 4
  plot(x = rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)], y = luad_pvals$hr[seq(1, length(luad_pvals$hr), 3)], ylim = c(0, 5), xlim = c(20, 500), pch = NA, ylab = 'HR', xlab = '#ClonalNeo used as threshold', las = 1)
  #segments(x0=as.numeric(rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)]), y0=luad_pvals$cil[seq(1, length(rownames(luad_pvals)), 3)],
  #         x1=as.numeric(rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)]), y1=luad_pvals$cih[seq(1, length(rownames(luad_pvals)),3)])
  points(x = rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)], y = luad_pvals$hr[seq(1, length(luad_pvals$hr), 3)], col = luad_pvals$col[seq(1, length(luad_pvals$col), 3)], pch = 16)
  points(x = rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)], y = rep(5, length(seq(1, length(rownames(luad_pvals)), 3))) , col = luad_pvals$ncol, pch =15)
  points(x = as.numeric(summary(tmp$ClonalNeo))[c(2,3,5)], y = c(5,5,5), pch = 15, col = 'gold')
  abline(h = 1, col = 'red', lty = 2)
  text(x = 100, y = 4.9, labels = 'number patients in high neo group', cex = 0.7)
  
  luad_pvals$ncol_s <- rbPal(30)[as.numeric(cut(luad_pvals$nu_s,breaks = 30))]
  luad_pvals$col_s <- ifelse(luad_pvals$p_s < 0.05, yes = 'red', no = 'black')
  luad_pvals$cih_s[which(luad_pvals$cih_s > 4)] <- 4
  plot(x = rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)], y = luad_pvals$hr_s[seq(1, length(luad_pvals$hr_s), 3)], ylim = c(0, 4), xlim = c(20, 500), pch = NA, ylab = 'HR', xlab = '#SubclonalNeo used as threshold', las = 1)
  #segments(x0=as.numeric(rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)]), y0=luad_pvals$cil_s[seq(1, length(rownames(luad_pvals)), 3)],
  #         x1=as.numeric(rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)]), y1=luad_pvals$cih_s[seq(1, length(rownames(luad_pvals)),3)])
  points(x = rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)], y = luad_pvals$hr_s[seq(1, length(luad_pvals$hr_s), 3)], col = luad_pvals$col_s[seq(1, length(luad_pvals$col_s), 3)], pch = 16)
  points(x = rownames(luad_pvals)[seq(1, length(rownames(luad_pvals)), 3)], y = rep(5, length(seq(1, length(rownames(luad_pvals)), 3))) , col = luad_pvals$ncol_s, pch =15)
  points(x = as.numeric(summary(tmp$SubclonalNeo))[c(2,3,5)], y = c(5,5,5), pch = 15, col = 'gold')
  abline(h = 1, col = 'red', lty = 2)
  text(x = 100, y = 4.9, labels = 'number patients in high neo group', cex = 0.7)
  
  luad_pvals <- keep_luad
  
  # # lusc
  tmp        <- patient.summary[which(patient.summary$simple.hist == 'LUSC'),]
  lusc_pvals <- data.frame(row.names = 1:max(tmp$ClonalNeo))
  lusc_pvals$hr <- NA
  lusc_pvals$p  <- NA
  lusc_pvals$cil <- NA
  lusc_pvals$cih <- NA
  lusc_pvals$nu <- NA
  lusc_pvals$hr_s <- NA
  lusc_pvals$p_s  <- NA
  lusc_pvals$cil_s <- NA
  lusc_pvals$cih_s <- NA
  lusc_pvals$nu_s <- NA
  for(n in 1:max(tmp$ClonalNeo)){
    test <- summary(coxph(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ ClonalNeo > n, data = tmp))
    hr  <- test$conf.int[1]
    cil <- test$conf.int[3]
    cih <- test$conf.int[4]
    p   <- test$sctest['pvalue']
    nu  <- survfit(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ ClonalNeo > n, data = tmp)$n[2]
    test <- summary(coxph(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ SubclonalNeo > n, data = tmp))
    hr_s <- test$conf.int[1]
    cil_s <- test$conf.int[3]
    cih_s <- test$conf.int[4]
    p_s  <- test$sctest['pvalue']
    nu_s <- survfit(Surv(allan_final_dsf_time, allan_final_dsf == "1") ~ SubclonalNeo > n, data = tmp)$n[2]
    lusc_pvals[as.character(n),'hr'] <- hr
    lusc_pvals[as.character(n),'p']  <- p
    lusc_pvals[as.character(n),'cil']  <- cil
    lusc_pvals[as.character(n),'cih']  <- cih
    lusc_pvals[as.character(n),'nu'] <- nu
    lusc_pvals[as.character(n),'hr_s'] <- hr_s
    lusc_pvals[as.character(n),'p_s']  <- p_s
    lusc_pvals[as.character(n),'cil_s']  <- cil_s
    lusc_pvals[as.character(n),'cih_s']  <- cih_s
    lusc_pvals[as.character(n),'nu_s'] <- nu_s
  }
  keep_lusc <- lusc_pvals
  lusc_pvals <- lusc_pvals[20:500,]
  
  
  lusc_pvals      <- keep_lusc[!is.na(keep_lusc$nu),]
  lusc_pvals$ncol <- rbPal(30)[as.numeric(cut(lusc_pvals$nu,breaks = 30))]
  lusc_pvals$col <- ifelse(lusc_pvals$p < 0.05, yes = 'red', no = 'black')
  lusc_pvals$cih[which(lusc_pvals$cih > 2)] <- 2
  plot(x = rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)], y = lusc_pvals$hr[seq(1, length(lusc_pvals$hr), 3)], ylim = c(0, 2), xlim = c(33, 474), pch = NA, ylab = 'HR', xlab = '#ClonalNeo used as threshold', las = 1)
  #segments(x0=as.numeric(rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)]), y0=lusc_pvals$cil[seq(1, length(rownames(lusc_pvals)), 3)],
  #         x1=as.numeric(rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)]), y1=lusc_pvals$cih[seq(1, length(rownames(lusc_pvals)),3)])
  points(x = rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)], y = lusc_pvals$hr[seq(1, length(lusc_pvals$hr), 3)], col = lusc_pvals$col[seq(1, length(lusc_pvals$col), 3)], pch = 16)
  points(x = rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)], y = rep(2.5, length(seq(1, length(rownames(lusc_pvals)), 3))) , col = lusc_pvals$ncol, pch =15)
  points(x = as.numeric(summary(tmp$ClonalNeo))[c(2,3,5)], y = c(2.5,2.5,2.5), pch = 15, col = 'gold')
  abline(h = 1, col = 'red', lty = 2)
  text(x = 120, y = 2.4, labels = 'number patients in high neo group', cex = 0.7)
  
  lusc_pvals      <- keep_lusc[!is.na(keep_lusc$nu_s),]
  lusc_pvals$ncol_s <- rbPal(30)[as.numeric(cut(lusc_pvals$nu_s,breaks = 30))]
  lusc_pvals$col_s <- ifelse(lusc_pvals$p_s < 0.05, yes = 'red', no = 'black')
  lusc_pvals$cih_s[which(lusc_pvals$cih_s > 5)] <- 5
  plot(x = rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)], y = lusc_pvals$hr_s[seq(1, length(lusc_pvals$hr_s), 3)], ylim = c(0, 8), xlim = c(0, 355), col = lusc_pvals$col_s[seq(1, length(lusc_pvals$col_s), 3)], pch = 16, ylab = 'HR', xlab = '#SubclonalNeo used as threshold', las = 1)
  #segments(x0=as.numeric(rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)]), y0=lusc_pvals$cil_s[seq(1, length(rownames(lusc_pvals)), 3)],
  #         x1=as.numeric(rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)]), y1=lusc_pvals$cih_s[seq(1, length(rownames(lusc_pvals)),3)])
  points(x = rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)], y = lusc_pvals$hr_s[seq(1, length(lusc_pvals$hr_s), 3)], col = lusc_pvals$col_s[seq(1, length(lusc_pvals$col_s), 3)], pch = 16)
  points(x = rownames(lusc_pvals)[seq(1, length(rownames(lusc_pvals)), 3)], y = rep(6, length(seq(1, length(rownames(lusc_pvals)), 3))) , col = lusc_pvals$ncol_s, pch =15)
  points(x = as.numeric(summary(tmp$SubclonalNeo))[c(2,3,5)], y = c(6,6,6), pch = 15, col = 'gold')
  abline(h = 1, col = 'red', lty = 2)
  text(x = 100, y = 5.9, labels = 'number patients in high neo group', cex = 0.7) 
}

# Figure S9F
{
  
  tmp <- patient.summary[which(patient.summary$simple.hist %in% c('LUAD', 'LUSC')),]
  tmp$low_evasion_clonal <- ifelse(tmp$low_evasion == TRUE | tmp$high_clonal == TRUE, yes =TRUE, no = FALSE)
  tmp$histology <- tmp$simple.hist
  
  tmp <- within(tmp, {
    histology <- factor(histology)
    sex <- factor(sex, labels = c("Female", "Male"))
    stage <- factor(allan_stage_no_3b, levels = c('1a', '1b', '2a', '2b', '3a'), labels = c('1a', '1b', '2a', '2b', '3'))
    Adjuvant.therapy <- factor(Adjuvant.therapy, levels = c('No adjuvant treatment', 'Adjuvant'))
    low_evasion <- factor(low_evasion, levels = c('FALSE', 'TRUE'))
    low_evasion_clonal <- factor(low_evasion_clonal, levels = c('FALSE', 'TRUE'))
    number.unique.hlas <- as.numeric(number.unique.hlas)
  })
  
  tmp_model <- coxph(Surv(allan_final_dsf_time, allan_final_dsf == '1') ~ histology + age + sex + stage + packyears + Adjuvant.therapy + low_evasion_clonal, data=tmp)
  ggforest(tmp_model, data = NULL, main = "Hazard ratio", fontsize = 0.7, refLabel = "reference", noDigits = 2)
}



tmp_model <- coxph(Surv(allan_final_dsf_time, allan_final_dsf == '1') ~ histology + age + sex + stage + packyears + Adjuvant.therapy + low_evasion + number.unique.hlas, data=tmp)
ggforest(tmp_model, data = NULL, main = "Hazard ratio", fontsize = 0.7, refLabel = "reference", noDigits = 2)
