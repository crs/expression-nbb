#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-e", "--expression"), action="store", type="character", help="Preprocessed expression file"),
  make_option(c("-s", "--sampleinfo"), action="store", type="character", help="Sample Info file"),
  make_option(c("-o", "--output"), action="store", type="character", help="Unified sample info output file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$output) | is.null(opt$sampleinfo) | is.null(opt$expression)) {
  stop("Expression, sampleinfo, and output filenames must be provided. See script usage (--help) for further information")
}


expr <- read.table(opt$expression, header=T, sep="\t")

expression.samples <- read.table(opt$sampleinfo, header=T)
expression.samples <- subset(expression.samples, Trivia_Name %in% colnames(expr))
expression.samples <- expression.samples[order(expression.samples$Trivia_Name),]

pmd <- sub('\\d*_(\\d*)', '\\1', expression.samples$PMD.h)
expression.samples$pmd  <- as.integer(pmd)


samplesheet <- data.frame(gender=expression.samples$Gender02, diagnosis=toupper(expression.samples$Diagnosis02), braak=expression.samples$Braak_Stage02, age=expression.samples$Age, apoe=expression.samples$ApoE, pmd=expression.samples$pmd)
rownames(samplesheet) <- expression.samples$Trivia_Name

write.table(samplesheet, quote=F,file=opt$output, row.names=T,col.names=T, sep=",")


