


rule all:
  input : 'stats/MUC.Stats.txt', 
          expand('plots/{cohort}-{gene}-{covar}.png', cohort=['MUC','NBB'], gene=['LASS1','GDF1'], covar=['diagnosis','braak']), 
          expand('plots/{cohort}-{gene}-{covar}.png', cohort=['MUC','NBB'], gene=['ILMN_1673232', 'ILMN_1728642', 'ILMN_1662280','ILMN_1685319'], covar=['diagnosis','braak']), 
          #'processed/MUC.neqc.filtered.txt',
          'processed/NBB.covariates.txt', 'stats/NBB.Covar.Correlations.txt',
          #expand('plots/{cohort}-{A}x{B}-{covar}.png',A=['ILMN_1728642'],B=['ILMN_1673232'], cohort=['MUC','NBB'], covar=['diagnosis', 'braak'])#,'stats/NBB.Stats.txt', 


rule plotRatio:
  input: expression='processed/{cohort}.neqc.filtered.transcript.txt', sampleinfo='processed/{cohort}.neqc.filtered.sampleinfo.txt'
  output: expand('plots/{{cohort}}-{{A}}x{{B}}-{covar}.png',covar=['diagnosis','braak'])  
  shell: './plotRatio.R -e {input.expression} -m {wildcards.A} -n {wildcards.B} -s {input.sampleinfo} -c diagnosis,braak -o plots/{wildcards.cohort}'

rule createPlots2:
  input: expression='processed/{cohort}.neqc.filtered.transcript.txt', sampleinfo='processed/{cohort}.neqc.filtered.sampleinfo.txt'
  output: expand('plots/{{cohort}}-{gene}-{covar}.png', gene=['ILMN_1673232','ILMN_1728642', 'ILMN_1662280','ILMN_1685319'], covar=['diagnosis','braak'])
  shell: './plotMarkers.R -e {input.expression} -m ILMN_1673232,ILMN_1728642,ILMN_1662280,ILMN_1685319 -s {input.sampleinfo} -c diagnosis,braak -o plots/{wildcards.cohort}'
 

rule createPlots:
  input: expression='processed/{cohort}.neqc.filtered.txt', sampleinfo='processed/{cohort}.neqc.filtered.sampleinfo.txt'
  output: expand('plots/{{cohort}}-{gene}-{covar}.png', gene=['LASS1','GDF1'], covar=['diagnosis','braak'])
  shell: './plotMarkers.R -e {input.expression} -m LASS1,GDF1 -s {input.sampleinfo} -c diagnosis,braak -o plots/{wildcards.cohort}'
  
rule correlateGenes:
  input: expr='processed/NBB.neqc.filtered.txt', covar='processed/NBB.covariates.txt'
  output: covarcor='stats/NBB.Covar.Correlations.txt', minmax='stats/NBB.Covar.Minmax.txt', genecorr='stats/NBB.Gene.Corr.txt'
  shell: './correlateGenes.R -e {input.expr} -c {input.covar} -o {output.covarcor} -m {output.minmax} -g {output.genecorr}'

rule alignCovars:
  input: uni='processed/NBB.neqc.filtered.sampleinfo.txt', covars='NBB/NBB.GWAS.Covariates.txt', samplesheet='NBB/NBB.GWAS.SamplesTable.txt'
  output: 'processed/NBB.covariates.txt'
  shell: './alignCovarsWithSamplesheet.R -s {input.uni} -c {input.covars} -g {input.samplesheet} -o {output}'

rule getStats:
  input: 'processed/{cohort}.neqc.filtered.sampleinfo.txt'
  output: 'stats/{cohort}.Stats.txt'
  shell: './getStats.R -s {input} -o {output}'

rule getMUCSampleInfo:
  input: expression='processed/MUC.neqc.filtered.txt', sampleinfo='MUC/MUC.SampleInfo.txt'
  output: 'processed/MUC.neqc.filtered.sampleinfo.txt'
  shell: './getMUCSampleInfoForPaper.R -e {input.expression} -s {input.sampleinfo} -o {output}'  

rule getNBBSampleInfo:
  input: expression='processed/NBB.neqc.filtered.txt', sampleinfo='NBB/NBB.SampleInfo.txt', samplesheet='NBB/NBB.Samplesheet.csv'
  output: 'processed/NBB.neqc.filtered.sampleinfo.txt'
  shell: './getNBBSampleInfoForPaper.R -e {input.expression} -s {input.sampleinfo} -t {input.samplesheet} -o {output}'

rule processMUC:
  input: probes='MUC/MUC.Sample.Probe.txt', control="MUC/MUC.Control.txt", dups="MUC/MUC.Duplicates.txt", rem="MUC/MUC.Remove.txt"
  output: avereps='processed/MUC.neqc.filtered.txt', transcript='processed/MUC.neqc.filtered.transcript.txt'
  shell: "./preprocessExpressionData.R -e {input.probes} -c {input.control} -d {input.dups} -r {input.rem} -o {output.avereps} -t {output.transcript}"
  
rule processNBB:
  input: probes='NBB/NBB.Sample.Probe.txt', dups="NBB/NBB.Duplicates.txt", rem="NBB/NBB.Remove.txt"
  output: avereps='processed/NBB.neqc.filtered.txt', transcript='processed/NBB.neqc.filtered.transcript.txt'
  shell: "./preprocessExpressionData.R -e {input.probes} -d {input.dups} -r {input.rem} -o {output.avereps} -x 1 -y 8 -t {output.transcript}"