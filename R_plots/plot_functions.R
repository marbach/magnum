# =============================================================================
# plot_functions.R
#
# Functions used by example_script.R
#
# The script is provided "as is" under the MIT license included with the magnum
# distribution.
#
# --
# Daniel Marbach
# April 27, 2015
# =============================================================================

require("ggplot2")
require("grid")
require("stringr")
require("data.table")


# =============================================================================
# (1) Load AUCs and compute p-values / scores

# -----------------------------------------------------------------------------
# Load AUC files for all given gwas and kernels in a matrix
loadAllAUCs <- function(dir, gwas, kernels, numPermut, k, logScale) {
    
    AUCs <- vector("list", length(kernels))
    #i<-1
    for (i in 1:length(kernels)) {
        kernel_i <- kernels[i]
        cat(i, "\t", kernel_i, "\n")
        
        # Matrix for all GWAS for kernel_i
        auc <- data.table(matrix(nrow=numPermut+1, ncol=length(gwas)))
        setnames(auc, gwas)
        
        #g <- 1
        for (g in 1:length(gwas)) {
            gwas_i <- gwas[g]
            
            # Load this specific gwas-kernel result
            auc_i <- loadAUC(dir, gwas_i, kernel_i, k, logScale)
            auc[, g:=auc_i, with=FALSE]
        }
        AUCs[[i]] <- auc
    }
    names(AUCs) <- kernels
    return(AUCs)
}


# -----------------------------------------------------------------------------
# Load a single AUC file
loadAUC <- function(dir, gwas_i, kernel_i, k, logScale) {
    
    filename <- paste(dir, "/", gwas_i, "--", kernel_i, ".AUC", sep="")
    filename.gz <- paste(filename, ".txt.gz", sep="")
    filename.txt <- paste(filename, ".txt", sep="")
    if (!file.exists(filename.gz)) {
        cat("File not found:", filename.gz, "\n")
        return(NA)
    }
    
    try(system(paste("gunzip -k", filename.gz)))
    #aucFile <- gzfile(filename)
    #auc <- read.table(aucFile, header=TRUE, colClasses="numeric")
    auc <- fread(filename.txt)
    file.remove(filename.txt)
    
    if (k == 5) {
        if (logScale)
            return(auc[, 10, with=FALSE])
        else
            return(auc[, 9, with=FALSE])
    }
    if (logScale)
        return(auc[, 4+k, with=FALSE])
    else
        return(auc[, k, with=FALSE])
}


# -------------------------------------------------------------------------------
# Compute empirical pvals based on AUCs
computeEnrichmentPvals <- function(aucs) {
    #aucs <- aucs.global
    pval = NULL
    for (i in 1:length(aucs)) {
        cat(i, "\n")
        #i <- 1
        auc_i = aucs[[i]]
        pval_i = NULL
        
        numSamples <- nrow(auc_i) - 1
        for (j in 1:ncol(auc_i)) {
            #j <- 1
            numGreater <- sum(auc_i[-1, j, with=FALSE] >= as.numeric(auc_i[1, j, with=FALSE]))
            pval_i[j] = (numGreater + 1) / (numSamples + 1)
        }
        pval <- rbind(pval, pval_i)
    }
    rownames(pval) <- names(aucs)  
    colnames(pval) <- colnames(aucs[[1]])
    return(pval)
}


# -------------------------------------------------------------------------------
# Correct pvals
fdrCorrection <- function(pvals) {
    
    # FDR correction
    p <- as.vector(pvals)
    p <- p.adjust(p, method="BH")
    
    pvals[,] <- p
    return(pvals)
}


# =============================================================================
#  (2) Plot connectivity enrichment curves

# -----------------------------------------------------------------------------
# Generate plots for all networks and parameters
generateEnrichmentCurvePlots <- function(dir, gwas, kernels,
                                         xrange=NULL, yrange=NULL,
                                         curveCutoff, pvals,
                                         outfile) {
  
  # Start pdf device
  numCol <- length(kernels)
  numRow <- length(gwas)
  cat("Writing file: ", outfile, "\n")
  pdf(file=paste(outfile, sep=""), height=numRow*8/2, width=numCol*7/2)  
  
  # Setup ggplot2
  theme_set(theme_classic())
  # Setup multi plot with grid
  vp.setup(numRow, numCol)
  
  significantGenes <- data.frame(gwas=NA, kernel=NA, qval=NA, numSignificantScores=NA)
  
  #g <- 1
  #k <- 1
  for (g in 1:length(gwas)) {
    for (k in 1:length(kernels)) {
      gwas_i <- gwas[g]
      kernel_i <- kernels[k]
      pval_i <- pvals[kernel_i, gwas_i]
      
      result <- enrichmentCurvePlot(dir, gwas_i, kernel_i, g, k,
                                    xrange, yrange,
                                    curveCutoff,
                                    pval_i)
      
      significantGenes <- rbind(significantGenes, result)
    }
  }
  # Turn off device, save
  dev.off()
  
  significantGenes <- significantGenes[-1,]  
  return(significantGenes)
}


# -----------------------------------------------------------------------------
# Plot enrichment curves for the given network and parameter
enrichmentCurvePlot <- function(dir, gwas_i, kernel_i, g, k,
                                xrange=NULL, yrange=NULL,
                                curveCutoff,
                                pval_i,
                                outfile) {
  #gwas_i <- gwas[1]
  #kernel_i <- kernels[1]
  #pval_i <- pvals[kernel_i, gwas_i]
  
  # Read the file
  filename <- paste(dir, "/", gwas_i, "--", kernel_i, sep="")
  curveFile <- gzfile(paste(filename, ".txt.gz", sep=""))
  curves <- read.table(curveFile, header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE, check.names=FALSE, colClasses="numeric")
  curvesEnrich <- curves[,1:9]
  
  # Compute fold change relative to expected
  curvesEnrich[,-c(1,4,9)] <- log2(curvesEnrich[,-c(1,4,9)] / curvesEnrich$median)
    
  # The genome-wide significance level for genes (Bonferroni threshold)
  numGenes <- (1/curveCutoff) * curves$k[nrow(curves)]
  genomeWideSignificanceLevel <- 0.05 / numGenes
  
  #g <- 1
  #k <- 1
  #title <- paste(getShortStudyNames(gwas_i))#, " (q=", round(10000*pval_i)/10000, ")", sep="")
  title <- gwas_i

  index <- g
  enrichmentCurvePlot.inner(curves=curvesEnrich, gwas_i, kernel_i, plotRow=index, plotCol=k,
                            xrange, yrange,
                            genomeWideSignificanceLevel, pval_i,
                            title)
  
  
  # Return the number of significant genes (scores, network clustering)
  isSignificant <- which(curves$geneScore <= genomeWideSignificanceLevel)
  numSignificantScores <- curves$k[isSignificant[length(isSignificant)]]  
  if (length(numSignificantScores) == 0)
      numSignificantScores = 0
  
  result <- data.frame(gwas=gwas_i, kernel=kernel_i, qval=pval_i, numSignificantScores=numSignificantScores)
  return(result)
}


# -----------------------------------------------------------------------------
# Plot enrichment curves (standard or window)
enrichmentCurvePlot.inner <- function(curves, gwas_i, kernel_i, plotRow, plotCol,
                                      xrange=NULL, yrange=NULL,
                                      genomeWideSignificanceLevel, pval_i,
                                      title) {
  
  # Initialize subplot
  subplot <- ggplot() + xlab("Ranked genes") + ylab("Enrichment")
  if (!is.null(title))
    subplot <- subplot + ggtitle(title)
  
  # Find the x-range
  obsCurves <- "observed"
  
  # x scale
  if (is.null(xrange))
    xrange <- range(curves$k)
  subplot <- subplot + scale_x_continuous(limits=xrange, expand=c(0,0))
  
  # y scale
  if (!is.null(yrange)) {
    subplot <- subplot + scale_y_continuous(limits=yrange, expand=c(0,0))
    
    # Truncate curves at ymin, ymax
    curves$observed <- truncateMe(curves$observed, yrange[1], yrange[2])
    for (i in 5:8)
      curves[,i] <- truncateMe(curves[,i], yrange[1], yrange[2])
    
  } else {
    subplot <- subplot + scale_y_continuous(expand=c(0,0))
  }
  
  # Plot fdr curves
  subplot <- subplot + 
      geom_ribbon(data=curves, mapping=aes(x=k, ymin=p_0.01, ymax=p_0.05), fill=rgb(171/255, 176/255, 211/255)) +
      geom_ribbon(data=curves, mapping=aes(x=k, ymin=p_0.05, ymax=p_0.95), fill=rgb(128/255, 134/255, 193/255)) +
      geom_ribbon(data=curves, mapping=aes(x=k, ymin=p_0.95, ymax=p_0.99), fill=rgb(171/255, 176/255, 211/255))
  
  # Plot observed enrichment
  data <- data.frame(x=curves$k, y=curves[,obsCurves])    
    
  subplot <- subplot + geom_line(data=data, mapping=aes(x=x, y=y), size=1, colour="red")
  
  # Plot y=0
  subplot <- subplot + geom_hline(yintercept=0, colour=rgb(0.15,0.15,0.15))
  
  # Plot x=genomeWideSignificanceLevel
  isSignificant <- which(curves$geneScore <= genomeWideSignificanceLevel)
  numSignificantScores <- curves$k[isSignificant[length(isSignificant)]]
  #alpha <- curves$k[which(curves$geneScore > genomeWideSignificanceLevel)[1]]
  subplot <- subplot + geom_vline(xintercept=numSignificantScores, colour=rgb(0.01,0.01,0.01))
  
  # Print the subplot
  #print(subplot)
  print(subplot, vp=vp.layout(plotRow, plotCol))
}


# -----------------------------------------------------------------------------
# Plot enrichment curves (standard or window)
truncateMe <- function(x, min, max) {
  
  x[x < min] <- min
  x[x > max] <- max
  return(x)
}


# -------------------------------------------------------------------------------
# Create multi-plot setup (nrow, ncol)
vp.setup <- function(x,y){ 
    # create a new layout with grid 
    grid.newpage() 
    # define viewports and assign it to grid layout 
    pushViewport(viewport(layout=grid.layout(x,y)))
}


# -------------------------------------------------------------------------------
# Access subplot (row, col)
vp.layout <- function(x,y){ 
    viewport(layout.pos.row=x, layout.pos.col=y)
}

