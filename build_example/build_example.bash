#git clone https://github.com/h3abionet/h3agwas-examples

#ls h3agwas-examples/data/imputed/vcf/*.vcf.gz > listvcf
#Cmt=1
#for file in h3agwas-examples/data/imputed/vcf/*.vcf.gz
#do
#if [ $Cmt -eq 1 ]
#then
#zcat $file > all.vcf
#Cmt=2
#else 
#zcat $file |grep -v "#" >> all.vcf
#fi
#done
##bcftools merge -l listvcf  -0 -Oz -o all_merged.vcf.gz
#cat all.vcf|bgzip -c > all.vcf.gz

filelistvcf=listvcf
ls h3agwas-examples/data/imputed/vcf/*.vcf.gz > $filelistvcf

filepheno=h3agwas-examples/data/pheno/pheno_test.all

filevcf=h3agwas-examples/data/imputed/all.vcf.gz

filelistvcf=listvcf
ls h3agwas-examples/data/imputed/vcf/*.vcf.gz > $filelistvcf
filevcf=h3agwas-examples/data/imputed/all.vcf.gz
filepheno=h3agwas-examples/data/pheno/pheno_test.all
pheno=pheno_qt1,pheno_qt2
fasta=~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz
##
file_snp_rel=list_pos_relatdness
awk '{print $1"\t"$4"\t"$4"\t"$1":"$4}' h3agwas-examples/data/array_plk/array.bim > $file_snp_rel



git clone https://github.com/h3abionet/h3agwas-examples

filelistvcf=listvcf
ls h3agwas-examples/data/imputed/vcf/*.vcf.gz > $filelistvcf
filevcf=h3agwas-examples/data/imputed/all.vcf.gz
FilePheno=h3agwas-examples/data/pheno/pheno_test.all
pheno=pheno_qt1,pheno_qt2
fasta=~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz

file_snp_rel=list_pos_relatdness
awk '{print $1"\t"$4"\t"$4"\t"$1":"$4}' h3agwas-examples/data/array_plk/array.bim > $file_snp_rel


## no loco, format 1 vcf in plink, no dosage used, independant sample relatdness limited position to file_snp_rel
#~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco0_dosage0_vcfi --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir sr1_loco0_dosage0_vcfi --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0 

## loco, format 1 vcf in plink, no dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco1_dosage0_vcfi --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir sr1_loco1_dosage0_vcfi --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0
exit


## no loco, format multi vcf in plink, no dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco0_dosage0_vcfmulti --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco o   --output_dir sr1_loco1_dosage0_vcfmulti --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0

## loco, format multi vcf in plink, no dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco0_dosage0_vcfmulti --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir sr1_loco1_dosage0_vcfi --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0

## no loco, format 1  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco0_dosage1_vcfi --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir sr1_loco0_dosage1_vcfi --file_vcf $filevcf  --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1

## loco, format 1  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco1_dosage1_vcfi --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir sr1_loco1_dosage1_vcfi --file_vcf $filevcf  $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1


## no loco, format multi  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco0_dosage1_vcfmulti --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir sr1_loco0_dosage1_vcfmulti --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1

## loco, format multi  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output sr1_loco1_dosage1_vcfmulti --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir sr1_loco1_dosage1_vcfmulti --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1


