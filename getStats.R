#!/usr/bin/Rscript

library(psych)
library(optparse)


option_list <- list(
  make_option(c("-s", "--sampleinfo"), action="store", type="character", help="Unified Sample Info file"),
  make_option(c("-o", "--output"), action="store", type="character", help="Output table")
)

opt <- parse_args(OptionParser(option_list=option_list))

expression.samples <- read.table(opt$sampleinfo, header=T, sep=',')

cases <- subset(expression.samples, diagnosis == 'AD')
controls <- subset(expression.samples, diagnosis == 'CTRL')


getApoe4Carrier <- function(carriers) {
  n.total <- length(carriers)
  n.carrier <- length(which(grepl('4', carriers)))
  return (n.carrier / n.total) * 100
}


getStats <- function(cohort) {
  age <- describe(cohort$age)
  pmd <- describe(cohort$pmd)
  braak <- describe(cohort$braak)
  apoe4.fraction <- getApoe4Carrier(cohort$apoe)
  
  return(data.frame(Age.Mean=age$mean, Age.SD=age$sd, PMD.mean=pmd$mean, PMD.sd=pmd$sd, Braak.Mean=braak$mean, Braak.SD=braak$sd, ApoE.Fraction=apoe4.fraction)) 
}


stats <- data.frame(samples=character(), N=integer(),Age.Mean=double(), Age.SD=double(), PMD.mean=double(), PMD.sd=double(), Braak.Mean=double(), Braak.SD=double(), ApoE.Fraction=double())

all.stats <- getStats(expression.samples)
all.stats <- cbind(data.frame(samples='All', N=nrow(expression.samples)), all.stats)
stats <- rbind(stats, all.stats)

print(stats)

cases.stats <- getStats(cases)
print(cases.stats)
stats <- rbind(stats, cbind(data.frame(samples='AD', N=nrow(cases)), cases.stats))
print(stats)

cases.male <- subset(cases, gender == 'M')
cases.male.stats <- getStats(cases.male)
stats <- rbind(stats, cbind(data.frame(samples='AD male', N=nrow(cases.male)), cases.male.stats))

cases.female <- subset(cases, gender == 'F')
cases.female.stats <- getStats(cases.female)
stats <- rbind(stats, cbind(data.frame(samples='AD female', N=nrow(cases.female)), cases.female.stats))

control.stats <- getStats(controls)
stats <- rbind(stats, cbind(data.frame(samples='Ctrl', N=nrow(controls)), control.stats))

control.male <- subset(controls, gender == 'M')
control.male.stats <- getStats(control.male)
stats <- rbind(stats, cbind(data.frame(samples='Ctrl male', N=nrow(control.male)), control.male.stats))

control.female <- subset(controls, gender == 'F')
control.female.stats <- getStats(control.female)
stats <- rbind(stats, cbind(data.frame(samples='Ctrl female', N=nrow(control.female)), control.female.stats))

write.table(stats, quote=F,file=opt$output, row.names=F,col.names=T, sep=",")