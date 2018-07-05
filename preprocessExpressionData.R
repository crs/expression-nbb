#!/usr/bin/Rscript


args <- commandArgs(trailingOnly = T)


## check correct script call
if (length(args) == 1000) {
  allArgs <- commandArgs()  
  pat <- "--file=./(.*)"
  scriptName <- sub(pat, "\\1", allArgs[grep("--file", allArgs)])
  cat("Usage: ", scriptName, " <expressionfile> (<controls>)\n")
  quit("no", 1, FALSE)
}

## 
library(optparse)
library(limma)

## parse options

option_list <- list(
  make_option(c("-e", "--expression"), action="store", type="character", help="Illumina expression file"),
  make_option(c("-c", "--control"), action="store", type="character", help="Illumina control file (optional)"),
  make_option(c("-d", "--duplicates"), action="store", type="character", help="File that contains duplicate identifiers (optional)"),
  make_option(c("-p", "--detection.p"), action="store", type="double", default=0.05, help="Detection p-value (default: 0.05)"),
  make_option(c("-s", "--detection.samples"), action="store", type="integer", default=3, help="Minimal amount of samples to detect (default: 3)"),
  make_option(c("-r", "--remove"), action="store", type="character", help="File that contains identifiers to remove (optional)"),
  make_option(c("-o", "--output"), action="store", type="character", help="Processed expression output file"),
  make_option(c("-t", "--transcriptoutput"), action="store", type="character", help="Processed transcript output file"),
  make_option(c("-x", "--sample_substring_start"), action="store", type="integer", help="Sample substring start (optional)"),
  make_option(c("-y", "--sample_substring_stop"), action="store", type="integer", help="Sample substring stop (optional)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# print(opt)

if (is.null(opt$expression)) {
  stop("Expression file must be provided. See script usage (--help) for further information")
}

if (!is.null(opt$control)) {
  expr <- read.ilmn(files=opt$expression, ctrlfiles=opt$control, sep="\t") 
} else {
  expr <- read.ilmn(files=opt$expression)
}

expressed.x <- rowSums(expr$other$Detection < opt$detection.p) >= opt$detection.samples
expressed.y <- colSums(expr$other$Detection < opt$detection.p) >= opt$detection.samples


expr <- expr[expressed.x, expressed.y]

expr <- neqc(expr)
#E <- expr$E
E <- avereps(expr, ID=expr$genes$SYMBOL)
E <- E$E
TE <- expr$E

if (!is.null(opt$remove)) {
  removeIDs <- readLines(opt$remove)
  print(colnames(E))
  print(removeIDs)
  removeIdx <- which(colnames(E) %in% removeIDs)
  E <- E[,-removeIdx]
  TE <- TE[,-removeIdx]  
}

if (!is.null(opt$duplicates)) {
  dups <- readLines(opt$duplicates)  
  for (i in 1:length(dups)) {
    id <- dups[i]
    duplicates <- grep(id, colnames(E))
    #print(id)
    #print(colnames(E))
    cat("Found duplicates:\n")
    if (length(duplicates) > 1) {
      print(duplicates)
      dup.E <- E[,grep(id, colnames(E))]
      dup.TE <- TE[,grep(id, colnames(TE))]
      
      te.means <- data.frame(rowMeans(dup.TE))
      means <- data.frame(rowMeans(dup.E))
      
      colnames(means)[1] <- id
      rownames(means) <- rownames(E)
      
      colnames(te.means)[1] <- id
      rownames(te.means) <- rownames(TE)
      
      cat("HERE_")
      #print(means)
      E <- cbind(E, means)
      E <- E[,-duplicates]
      
      TE <- cbind(TE, te.means)
      TE <- TE[,-duplicates]
    }
  }
}

if (!is.null(opt$sample_substring_start) & !is.null(opt$sample_substring_stop)) {
  colnames(E) <- substr(colnames(E),opt$sample_substring_start, opt$sample_substring_stop)
  colnames(TE) <- substr(colnames(TE),opt$sample_substring_start, opt$sample_substring_stop)
}

print(E[1:10,1:10])
print(TE[1:10,1:10])

E <- E[,order(colnames(E))]
TE <- TE[,order(colnames(TE))]

### HACK 
#if (grepl('\\d*-\\d*', colnames(E))  
#colnames(E) <- paste('NBB-',colnames(E),sep="")

write.table(E, file=opt$output, quote=F, row.names=T, col.names=T, sep="\t")
write.table(TE, file=opt$transcriptoutput,  quote=F, row.names=T, col.names=T, sep="\t")
