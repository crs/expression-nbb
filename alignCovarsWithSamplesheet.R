#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-s", "--samplesheet"), action="store", type="character", help="Unified sample sheet"),
  make_option(c("-c", "--covars"), action="store", type="character", help="Covariate file with FID and IID"),
  make_option(c("-g", "--gwas"), action="store", type="character", help="GWAS Samplesheet (with 'Index' and 'Sample ID' columns) and SampleName"),
  make_option(c("-o", "--output"), action="store", type="character", help="Output covar file")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$samplesheet) | is.null(opt$covars) | is.null(opt$gwas)) {
  stop("Wrong parameters. See script usage (--help) for further information")
}

samples.table <- read.delim(opt$gwas)
covariates <- read.delim(opt$covars, sep=" ", header=T)

samples.table <- samples.table[order(samples.table$Index),]
samples.cov <- merge(samples.table[,c("SampleName", "Index", "Sample.ID", 'Gender', "Diagnosis")], covariates, by.x = c("Index","Sample.ID"), by.y=c("FID","IID"))

samples.cov$SampleName <- gsub('_','-', samples.cov$SampleName)

expression.samples <- read.table(opt$samplesheet, header=T, sep=',', check.names=F)

samples <- samples.cov[which(samples.cov$SampleName %in% rownames(expression.samples)),]

head(samples)

write.table(samples, file=opt$output, quote=F, row.names=T, col.names=T, sep="\t")

