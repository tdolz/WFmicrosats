library("diveRsity") #this one fucked everything. 

divRNI <-read.delim("/Users//tdolan/Documents//R-Github//WFmicrosats/divRNI")
#divRNIrare<-read.delim("/Users//tdolan/Documents//R-Github//WFmicrosats/divRNI_rareified")

t <-diveRsity::divMigrate(infile=divRNI, outfile = NULL, boots = 1000, stat = "all", 
                          filter_threshold = 0.3, plot_network = TRUE, 
                          plot_col = "darkblue", para = FALSE)

s <-divMIG(infile=divRNI, outfile = NULL, boots = 1000, stat = "all", 
                          filter_threshold = 0.3, plot_network = TRUE, 
                          plot_col = "darkblue", para = FALSE)






#create own distance matrix from the MIGRATE outfile

#all individuals outfile. 

allmig <-matrix(c(NA,1,0.995,1.053,0.993,1.027,NA,0.959,1.047,0.955,0.989,0.987,NA,1.116,0.917,1.001,1.015,1.013,NA,1.023,0.985,0.959,0.912,0.937,NA), nrow=5,ncol=5, byrow=TRUE)

noadultsmig <-matrix(c())

scaledmig <-allmig/1.116



#attempt to insert your matrix into the graph
qgraph::qgraph(allmig, nodeNames = sapply(dat$indnms, "[", 1),
                 legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
               mar = c(2, 2, 5, 5), curve = 2.5, threshold=0, edge.width=1, fade=TRUE, colFactor=10)

qgraph::qgraph(scaledmig, nodeNames = sapply(dat$indnms, "[", 1),
               legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
               mar = c(2, 2, 5, 5), curve = 2.5, threshold=0, edge.width=1, fade=TRUE, colFactor=10)


dRelPlt
allmig

qgraph::qgraph(dRelPlt, nodeNames = NULL, #sapply(dat$indnms, "[", 1),
               legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
               mar = c(2, 2, 5, 5), curve = 2.5, threshold=0, edge.width=1, fade=TRUE)

dRel <- dMig/max(dMig, na.rm = TRUE)
dRel[is.nan(dRel)] <- NA
dRelPlt <- dRel
dRelPlt[dRelPlt < filter_threshold] <- 0

plot_col <-"darkblue"
filter_threshold <- 0.5
allrel <-allmig/max(allmig, na.rm =TRUE)
allrel[is.nan(allrel)] <-NA
allrelPlt <-allrel
allrelPlt[allrelPlt < filter_threshold] <-0

qgraph::qgraph(allrelPlt, nodeNames = NULL, #sapply(dat$indnms, "[", 1),
               legend = TRUE, edge.labels = TRUE, 
               mar = c(2, 2, 5, 5), curve = 2.5)


######investigate divMIGRATE  ##########

#divMIG <- function (infile = NULL, outfile = NULL, boots = 0, stat = "all", 
 #         filter_threshold = 0, plot_network = FALSE, plot_col = "darkblue", 
      #    para = FALSE) 
#{
outfile=NULL
infile <- divRNI
stat = "all"
filter_threshold = 0
plot_network =TRUE
plot_col = "darkblue"
para =FALSE
boots = 0

  dat <- rgp(infile)
  npops <- length(dat$genos)
  nloci <- length(dat$af)
  dat$af <- lapply(dat$af, function(x) {
    cs <- colSums(x)
    x[, cs == 0] <- NA
    return(x)
  })
  if (!is.null(outfile)) {
    dir.create(path = paste(getwd(), "/", outfile, "-[divMigrate]", 
                            "/", sep = ""))
    of <- paste(getwd(), "/", outfile, "-[divMigrate]", "/", 
                sep = "")
  }
  pw <- combn(npops, 2)
  hths <- lapply(dat$af, pwHt, pw = pw - 1)
  ht <- lapply(hths, "[[", "ht")
  hs <- lapply(hths, "[[", "hs")
  if (stat == "d" || stat == "all" || stat == "Nm") {
    d <- function(ht, hs) {
      return(((ht - hs)/(1 - hs)) * 2)
    }
  }
  if (stat == "gst" || stat == "all" || stat == "Nm") {
    g <- function(ht, hs) {
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
  }
  if (stat == "Nm" || stat == "all") {
    Nm <- function(g, d, n) {
      t1 <- (1 - g)/g
      t2 <- ((n - 1)/n)^2
      t3 <- ((1 - d)/(1 - ((n - 1)/n) * d))
      return(0.25 * t1 * t2 * t3)
    }
  }
  if (stat == "d" || stat == "all" || stat == "Nm") {
    dloc <- mapply(d, ht = ht, hs = hs, SIMPLIFY = "array")
    dloc[is.nan(dloc)] <- 1
    hrmD <- apply(dloc, c(1, 2), function(x) {
      mn <- mean(x, na.rm = TRUE)
      vr <- var(x, na.rm = TRUE)
      return(1/((1/mn) + vr * (1/mn)^3))
    })
    dMig <- (1 - hrmD)/hrmD
    dMig[is.infinite(dMig)] <- NA
    dRel <- dMig/max(dMig, na.rm = TRUE)
    dRel[is.nan(dRel)] <- NA
  }
  if (stat == "gst" || stat == "all" || stat == "Nm") {
    g <- function(ht, hs) {
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
    hsAr <- array(unlist(hs), dim = c(npops, npops, nloci))
    mnHs <- apply(hsAr, c(1, 2), mean, na.rm = TRUE)
    htAr <- array(unlist(ht), dim = c(npops, npops, nloci))
    mnHt <- apply(htAr, c(1, 2), mean, na.rm = TRUE)
    hrmGst <- g(mnHt, mnHs)
    gMig <- ((1/hrmGst) - 1)/4
    gMig[is.infinite(gMig)] <- NA
    gRel <- gMig/max(gMig, na.rm = TRUE)
  }
  if (stat == "all" || stat == "Nm") {
    nm <- Nm(hrmGst, hrmD, 2)
    diag(nm) <- NA
    nmRel <- nm/max(nm, na.rm = TRUE)
  }
  if (plot_network || boots != 0L) {
    if (stat == "d" || stat == "all" || stat == "Nm") {
      dRelPlt <- dRel
      dRelPlt[dRelPlt < filter_threshold] <- 0
    }
    if (stat == "gst" || stat == "all" || stat == "Nm") {
      gRelPlt <- gRel
      gRelPlt[gRelPlt < filter_threshold] <- 0
    }
    if (stat == "all" || stat == "Nm") {
      nmRelPlt <- nmRel
      nmRelPlt[nmRelPlt < filter_threshold] <- 0
    }
  }
  if (plot_network) {
    if (stat == "d") {
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; D)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; D)", sep = ""))
      }
    }
    if (stat == "gst") {
      qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Gst)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Gst)", sep = ""))
      }
    }
    if (stat == "Nm") {
      qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                  "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Nm)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Nm)", sep = ""))
      }
    }
    if (stat == "all") {
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; D)", sep = ""))
      qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Gst)", sep = ""))
      qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                  "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Nm)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; D)", sep = ""))
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Gst)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Nm)", sep = ""))
      }
    }
    if (boots == 0L && !is.null(outfile)) {
      dev.off()
    }
  }
  if (boots != 0L) {
    ps <- sapply(dat$indnms, length)
    idx <- lapply(1:boots, function(i) {
      lapply(ps, function(x) {
        return(sample(x, size = x, replace = TRUE))
      })
    })
    if (para) {
      cl <- parallel::makeCluster(detectCores())
      parallel::clusterExport(cl, c("bsFun", "dat", "pw", 
                                    "stat"), envir = environment())
      bsStat <- parallel::parLapply(cl, idx, function(x) {
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, 
                     pw = pw, stat = stat))
      })
      parallel::stopCluster(cl)
    }
    else {
      bsStat <- lapply(idx, function(x) {
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, 
                     pw = pw, stat = stat))
      })
    }
    if (stat == "d" || stat == "all") {
      bsD <- sapply(bsStat, "[[", "dRel", simplify = "array")
    }
    if (stat == "gst" || stat == "all") {
      bsG <- sapply(bsStat, "[[", "gRel", simplify = "array")
    }
    if (stat == "Nm" || stat == "all") {
      bsNm <- sapply(bsStat, "[[", "nmRel", simplify = "array")
    }
    sigDiff <- function(x, y) {
      if (x[1] < y[1] && x[2] < y[1]) {
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    }
    if (stat == "d" || stat == "all") {
      sigMatD <- matrix(NA, nrow = ncol(dRel), ncol(dRel))
      for (i in 1:ncol(pw)) {
        p1 <- quantile(bsD[pw[1, i], pw[2, i], ], prob = c(0.025, 
                                                           0.975))
        p2 <- quantile(bsD[pw[2, i], pw[1, i], ], prob = c(0.025, 
                                                           0.975))
        sigMatD[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
        sigMatD[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
      }
      dRelPlt[!sigMatD] <- 0
    }
    if (stat == "gst" || stat == "all") {
      sigMatG <- matrix(NA, nrow = ncol(gRel), ncol(gRel))
      for (i in 1:ncol(pw)) {
        p1 <- quantile(bsG[pw[1, i], pw[2, i], ], prob = c(0.025, 
                                                           0.975))
        p2 <- quantile(bsG[pw[2, i], pw[1, i], ], prob = c(0.025, 
                                                           0.975))
        sigMatG[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
        sigMatG[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
      }
      gRelPlt[!sigMatG] <- 0
    }
    if (stat == "Nm" || stat == "all") {
      sigMatNm <- matrix(NA, nrow = ncol(nmRel), ncol(nmRel))
      for (i in 1:ncol(pw)) {
        p1 <- quantile(bsNm[pw[1, i], pw[2, i], ], prob = c(0.025, 
                                                            0.975))
        p2 <- quantile(bsNm[pw[2, i], pw[1, i], ], prob = c(0.025, 
                                                            0.975))
        sigMatNm[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
        sigMatNm[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
      }
      nmRelPlt[!sigMatNm] <- 0
    }
    if (plot_network) {
      if (stat == "d" && !is.null(outfile)) {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
        dev.off()
      }
      if (stat == "gst" && !is.null(outfile)) {
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Gst method)", sep = ""))
        dev.off()
      }
      if (stat == "Nm" && !is.null(outfile)) {
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
        dev.off()
      }
      if (stat == "all" && !is.null(outfile)) {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
        dev.off()
      }
      if (stat == "d") {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
      }
      if (stat == "gst") {
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Gst method)", sep = ""))
      }
      if (stat == "Nm") {
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
      }
      if (stat == "all") {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Gst Method)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
      }
    }
  }
  if (boots != 0L) {
    if (stat == "d") {
      list(dRelMig = dRel, dRelMigSig = dRelPlt)
    }
    else if (stat == "gst") {
      list(gRelMig = gRel, gRelMigSig = gRelPlt)
    }
    else if (stat == "Nm") {
      list(nmRelMig = nmRel, nmRelMigSig = nmRelPlt)
    }
    else if (stat == "all") {
      list(dRelMig = dRel, dRelMigSig = dRelPlt, gRelMig = gRel, 
           gRelMigSig = gRelPlt, nmRelMig = nmRel, nmRelMigSig = nmRelPlt)
    }
  }
  else {
    if (stat == "d") {
      list(dRelMig = dRel)
    }
    else if (stat == "gst") {
      list(gRelMig = gRel)
    }
    else if (stat == "Nm") {
      list(nmRelMig = nmRel)
    }
    else if (stat == "all") {
      list(dRelMig = dRel, gRelMig = gRel, nmRelMig = nmRel)
    }
  }
}
<bytecode: 0x7f8b831d3ee0>
  <environment: namespace:diveRsity>