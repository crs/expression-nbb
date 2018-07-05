#!/usr/bin/Rscript

## check correct script call
if (length(args) == 1000) {
  allArgs <- commandArgs()  
  pat <- "--file=./(.*)"
  scriptName <- sub(pat, "\\1", allArgs[grep("--file", allArgs)])
  cat("Usage: ", scriptName, " <expressionfile> (<controls>)\n")
  quit("no", 1, FALSE)
}

cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}


## 
library(optparse)
library(limma)

## parse options

option_list <- list(
  make_option(c("-e", "--expression"), action="store", type="character", help="Processed expression file with marker in rows"),
  make_option(c("-m", "--markers"), action="store", type="character", help="Markers (comma separated)"),
  make_option(c("-c", "--strat"), action="store", type="character", help="Values to stratify against (comma separated)"),  
  make_option(c("-s", "--samples"), action="store", type="character", help="Sample Info"),
  make_option(c("-o", "--output"), action="store", type="character", help="Output prefix")
)

if (!exists('arguments')) {
  arguments = commandArgs(trailingOnly = T)
}

opt <- parse_args(OptionParser(option_list=option_list), args = arguments)

if (is.null(opt$expression) | is.null(opt$markers) | is.null(opt$output)) {
  stop("Expression file, markers and output prefix must be provided. See script usage (--help) for further information")
}

#print(opt$markers)
markers <- strsplit(opt$markers, ',')[[1]]
strats <- strsplit(opt$strat, ',')[[1]]

print(markers)

E <- read.table(opt$expression, header=T, sep="\t", check.names=F)
samples <- read.table(opt$samples,  header=T, sep=',')

#E[1:10,1:10]
#samples[1:5,1:5]

library(ade4)

for (marker in markers) {  
    for (strat in strats) {
      tryCatch({
        if (grepl('NBB', opt$expression)) {
          rin <- which(samples[,'rin'] == 1)
          gene.strat <- samples[rin,strat]
          gene.expr <- as.numeric(E[marker,rin])                    
        } else {
          gene.strat <- samples[,strat]
          gene.expr <- as.numeric(E[marker,])            
        }
        if (class(gene.strat) == 'numeric' | class(gene.strat) == 'integer') {
          gene.strat = as.factor(round(gene.strat))
        }
        
        #op <- par(mar = rep(0, 4))
        png(filename = paste(opt$output, '-', marker, '-', strat, '.png', sep=""), width = 900, height = 900, units = "px", pointsize = 12)
        #par(op)
		
        if (marker == 'ILMN_1673232' & strat == 'braak') {
            cases = which(as.integer(gene.strat) > 3)
            print(gene.strat[cases])
            gene.expr = gene.expr[cases]
            gene.strat = gene.strat[cases]
            plot(gene.expr ~ gene.strat, xlab=cap(strat), ylab=marker, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
            points(gene.expr ~ gene.strat)
        } else {        
            plot(gene.expr ~ gene.strat, xlab=cap(strat), ylab=marker, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
            points(gene.expr ~ gene.strat)
        }
        if (strat == 'diagnosis') {
          #print(gene.strat)
          #print(which(gene.strat == 'AD'))
          #print(which(gene.strat == 'CTRL'))
          test = t.test(gene.expr[which(gene.strat == 'AD')], gene.expr[which(gene.strat == 'CTRL')])
          #print(test)
          #print(test$p.value)
          #text(20, 20, labels = c('p-Val:', test$p.value), pos = 4)
          legend("bottomleft",legend=paste('p-val:', test$p.value), bty ="n", cex=2.0, pch=1, pt.cex = 1) 
        }
        if (strat == 'braak') {        
          model = lm(gene.expr ~ as.numeric(gene.strat))
          summ = summary(model)
          print(summ)
          print(names(summ))
          #print(names(summ$residuals))
          #pvalue = summ$coefficients[1,4]

          legend("bottomleft",legend=paste('adj R-squared', summ$adj.r.squared), bty ="n", cex=2.0, pch=1, pt.cex = 1) 
        }
        dev.off()
      }, 
      warning = function(w) {
          print(paste("Warning: Marker", marker, "not found"))
          png(filename = paste(opt$output, '-', marker, '-', strat, '.png', sep=""), width = 1200, height = 1200, units = "px", pointsize = 12)
          plot.new()
          text(paste("Warning: Marker", marker, "not found"))
          dev.off()
      }, 
      error = function(e) {
      }, 
      finally = {         
      })
  }
}
