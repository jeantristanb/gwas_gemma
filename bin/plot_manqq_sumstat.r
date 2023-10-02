#!/usr/bin/env Rscript
library("optparse")
library(data.table)
library(fastman)

option_list = list(
  make_option(c( "--data"), type="character",
              help="data summary statistics", metavar="character"),
  make_option(c("--p_header"), type="character",
              help="ped file contains genotype for each individual", metavar="character"),
  make_option(c( "--chr_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--bp_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--af_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--rs_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--maf"), 
              help="phenotype in data", type="numeric", default=0.01),
  make_option(c( "--out"), type="character",
              help="phenotype in data", metavar="character")
);

opt = OptionParser(option_list=option_list);
opt = parse_args(opt);
test=F
if(test){
opt=list("data"="gwas_merge/dbp/env12/dbp_all_env12.gwas","af_header"="EAF_ALL", "p_header"="P_INT_ROBUST","chr_header"="CHR",out= "gwas_merge/dbp/env12/dbp_all_env12")
}


alldata<-fread(opt[['data']])
out=opt[['out']]

png(paste(out,'_distchro.png',sep=''),   width = 480*3, height = 480*2, res=200)
plot(table(alldata[[opt[['chr_header']]]]), xlab='chromosome', ylab='distribution')
dev.off()
if(!is.null(opt[['af_header']])){
 maf=opt[['maf']];maf2<-1 - maf
 alldataqq<-alldata[alldata[[opt[['af_header']]]]>maf & alldata[[opt[['af_header']]]]<maf2,]
}else {
alldataqq<-alldata
}
listhead<-c(opt[['chr_header']], opt[['bp_header']], opt[['p_header']], opt[['rs_header']])
if(all(listhead %in% names(alldataqq))==F){
cat('header not found', paste(listhead[!(listhead %in% names(alldataqq))], collapse=','),'\nexit\n')
q('n')
}
alldataqq<-na.omit(alldataqq[,..listhead])

png(paste(out,'_qq.png',sep=''))
fastqq (alldataqq[[opt[['p_header']]]])
dev.off()

rsheader=opt[['rs_header']]
if(is.null(opt[['rs_header']])){
alldataqq$rsid<-paste(alldataqq[[opt[['chr_header']]]], alldataqq[[opt[['bp_header']]]], sep=':')
rsheader='rsid'
}else{
balise<-is.na(alldataqq[[rsheader]]) | alldataqq[[rsheader]]=='.' 
alldataqq[[rsheader]][balise]<-paste(alldataqq[[opt[['chr_header']]]][balise], alldataqq[[opt[['bp_header']]]][balise], sep=':')
}
print(alldataqq)
png(paste(out,'_man.png', sep=''),   width = 480*6, height = 480*2, res=200)
fastman(alldataqq,  chr =opt[['chr_header']] , bp = opt[['bp_header']], p = opt[['p_header']], snp=rsheader,annotatePval=5E-8)
dev.off()


