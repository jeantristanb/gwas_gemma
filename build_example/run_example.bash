scl_enabled java-11
head=$1

workdir="-w /spaces/jeantristan/gwas_ckgen/$head"


file_snp_rel=list_pos_relatdness
filelistvcf=listvcf
filevcf=h3agwas-examples/data/imputed/all.vcf.gz
FilePheno=h3agwas-examples/data/pheno/pheno_test_missing.all

pheno=pheno_qt1,pheno_qt1_miss
fasta=~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz


#head="sr1_loco0_dosage0_vcfi"

if [ $head == "sr1_loco0_dosage1_vcfi_covnores" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 --covariates Sex,batch --pheno_residuals 0  $workdir 
exit $?
fi

if [ $head == "sr1_loco0_dosage1_vcfi_covres" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 --covariates Sex,batch --pheno_residuals 1  $workdir 
exit $?
fi



## no loco, format 1 vcf in plink, no dosage used, independant sample relatdness limited position to file_snp_rel
if [ $head == "sr1_loco0_dosage0_vcfi" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0 $workdir
exit $?
fi

if [ $head == "sr1_loco0_dosage0_vcfi_addpc" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0  --addpcs 10 $workdir
exit $?
fi
#params.pheno_tr_fct=""
#params.phenores_tr_fct=""
#params.pheno_residuals=1

if [ $head == "sr1_loco0_dosage0_vcfi_addpc_nores" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0  --addpcs 10 --pheno_residuals 0 $workdir
exit $?
fi

if [ $head == "sr1_loco0_dosage0_vcfi_addpc_res_rin" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0  --addpcs 10 --pheno_residuals 1 --phenores_tr_fct invnorm $workdir
exit $?
fi

if [ $head == "sr1_loco0_dosage1_vcfi_addpc_res_rin" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1  --addpcs 10 --pheno_residuals 1 --phenores_tr_fct invnorm $workdir
exit $?
fi
## loco, format 1 vcf in plink, no dosage used, to sample relatdness limited position to file_snp_rel


if [ $head == "sr1_loco0_dosage0_vcfmulti" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco o   --output_dir $head --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0 $workdir
exit $?
fi

if [ $head == "sr1_loco0_dosage1_vcfi" ]
then
## no loco, format 1  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --file_vcf $filevcf  --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 $workdir
exit $?
fi

if [ $head == "sr1_loco1_dosage1_vcfi" ]
then
## loco, format 1  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir $head --file_vcf $filevcf  $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 $workdir
exit $?
fi


if [ $head == "sr1_loco0_dosage1_vcfmulti" ]
then 
## no loco, format multi  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 0   --output_dir $head --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 $workdir
exit $?
fi

## loco, format multi  vcf in plink/bimbam, dosage used, to sample relatdness limited position to file_snp_rel

if [  $head == "sr1_loco1_dosage1_vcfmulti" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir $head --listfile_vcf $filelistvcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 $workdir
exit $?
fi


if [  $head == "sr1_loco1_dosage0_vcfi" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir $head --file_vcf $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 0 $workdir
exit $?
fi


if [ ! -f h3agwas-examples/utils/listbimbam ]
then
rm -f h3agwas-examples/utils/listbimbam h3agwas-examples/utils/listbimbam_ind h3agwas-examples/utils/listbimbam_annot
for chro in `seq 1 22`
do
filebimbam=h3agwas-examples/data/imputed/bimbam_chro/$chro".pbwt"*.bimbam  
echo $chro" "$filebimbam >>h3agwas-examples/utils/listbimbam
bimbam_annot=h3agwas-examples/data/imputed/bimbam_chro/$chro".pbwt"*.annotation 
echo $chro" "$bimbam_annot >>h3agwas-examples/utils/listbimbam_annot
done
fi

if [  $head == "sr1_loco1_dosage1_bimbammulti" ]
then
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno $pheno  --output $head --gemma 1 -profile slurmSingularity  --sample_snps_rel 1  -resume --gemma_loco 1   --output_dir $head --file_ $filevcf --reffasta $fasta --snps_include_rel $file_snp_rel --dosage 1 $workdir --listfile_bimbam h3agwas-examples/utils/listbimbam  --file_bimbam_ind h3agwas-examples/data/imputed/bimbam_chro/1_sample_vcf.keep --bfile h3agwas-examples/data/imputed/imput_data.clean --gemma_mem_req_rel 20GB
exit $?
fi


