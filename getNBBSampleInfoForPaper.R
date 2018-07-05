#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-e", "--expression"), action="store", type="character", help="Preprocessed expression file"),
  make_option(c("-s", "--sampleinfo"), action="store", type="character", help="Sample Info file"),
  make_option(c("-o", "--output"), action="store", type="character", help="Unified sample info output file"),
  make_option(c("-t", "--samplesheet"), action="store", type="character", help="Samplesheet with RIN (optional)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$output) | is.null(opt$sampleinfo) | is.null(opt$expression)) {
  stop("Expression, sampleinfo, and output filenames must be provided. See script usage (--help) for further information")
}

expr <- read.table(opt$expression, header=T, sep="\t", check.names=F, nrows=1)

#head(expr,n=2)

expression.samples <- read.delim2(opt$sampleinfo, header=T, sep="\t", check.names=F)
expression.samples$NBB <- sub('_','-', expression.samples$NBB)
#expression.samples$NBB <- paste('NBB-', expression.samples$NBB, sep="")
expression.samples <- subset(expression.samples, NBB %in% colnames(expr))
expression.samples <- expression.samples[order(expression.samples$NBB),]

expression.samples$pmd  <- expression.samples$PMD_MIN / 60

if (opt$samplesheet != '') {
  samples_rin <- read.delim2(opt$samplesheet, skip=7, sep=',',stringsAsFactors=F)
  samples_rin$Sample_Name <- substr(samples_rin$Sample_Name,1,8)
  samples_rin <- subset(samples_rin, Sample_Name %in% expression.samples$NBB)
  samples_rin <- samples_rin[order(samples_rin$Sample_Name),]  
  samples_rin <- samples_rin[-which(duplicated(samples_rin$Sample_Name)),]    
  print(nrow(samples_rin))
  print(nrow(expression.samples))
}

samplesheet <- data.frame(gender=expression.samples$SEX, diagnosis=toupper(expression.samples$DIAGNOSIS), braak=expression.samples$BRAAK, age=expression.samples$AGE, apoe=expression.samples$APOE, pmd=expression.samples$pmd, rin=samples_rin$QualityScore)
rownames(samplesheet) <- expression.samples$NBB

write.table(samplesheet, quote=F,file=opt$output, row.names=T,col.names=T, sep=",")