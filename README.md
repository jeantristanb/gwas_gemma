# GWAS using GEMMA

## input :

### plink : 
 * `input_pat`
 * `input_dir`
  * `bfile` : 

### vcf :
 * `list_vcf` : ""
   * `file_listvcf`  : ""
 

### ref fasta

* `reffasta` [""]

### bim bam
 * bim bam dosage :
   * `file_bimbam` : do you have bim bam file
   * `listfile_bimbam` do you have a file ocntain list of bim bam
## formating vcf (option)
 * `min_scoreinfo`  [default : 0.3]
 * `genotype_field` : default : 
 * `score_imp` [default INFO]
 * `unzip_zip` : for vcf if there is zip and password [0]
 * `unzip_password` = ""


## filter 
 * `cut_maf`  [default : 0.01]
 * `keep_vcf` : keep individual on vcf 

## build relatdness
 * `listsnps_buildrelat` [default : ""]
 * Sample snps for relatdness :
  * `snps_include_rel` : snp that we used to defined snps of relatdness [ default : ""]
  * `snps_exclude_rel` :  snp that we used to exclude some snps of relatdness [default : ""] 
  * `sample_snps_rel`  :  do we need to smaple snps for reladness
  * `cut_maf_rel`  : [default : "0.01"]
  * `plink_indep_pairwise` [default : "100 20 0.1"]
  * `thin_snp_rel` : max snps for relatdnees ""

## bin
* `vcfftools_bin` [default : 'vcftools']
* `gemma_bin`  [default :"gemma"]
* `qctoolsv2_bin` ["qctool"]
* bcftools_bin [default "bcftools"]

## gemma parameter
 * performed gemma on dosage ? `dosage` [default: 0 ]
 *  `gemma_multi`  : split by chromosome
 * `gemma_mem_req` : memory for gemma "6GB"
 * `gemma_mem_req_rel` : memory for relatdness
 * `gemma_mat_rel` : [default : ""]
 * `gemma_num_cores` [cpus for gemma default : 8]
 * `gemma_loco`  : performed a loco with gemma yes : 1, no 0 [ default :0 ]

##  phenotype 
 * `pheno` : phenotype


## tocheck 
* `exclude_snps`
* `rs_list=`""

## other :
 * `other_process_mem_req`  [default '10G' ]
 * plink_mem_req = '6GB' // how much plink needs for this


# Example 
#nextflow pull h3abionet/h3agwas
## first test : gemma with plink no loco
```
git clone https://github.com/h3abionet/h3agwas-examples

ls h3agwas-examples//*.vcf.gz > listvcf
FilePheno=exampledata2b.pheno
awk '{$5=$4;print $0}' ~/Travail/git/h3agwas_buildexample/data2//exampledata2.pheno > exampledata2b.pheno

## no loco, sample relatdness 
~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_plk_vcf  --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --reffasta ~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz 

~/nextflow run gwas_ckdawigen/assoc.nf --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_plk_multivcf  --listfile_vcf listvcf --reffasta ~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz 

~/nextflow run gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_1vcf --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --dosage 1

zcat ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz | head -1000 |grep "#"  | tail -1 |awk '{for(Cmt=40;Cmt<NF;Cmt++)print $Cmt}' > filekeep_vcf
~/nextflow run gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_1vcf_indkeep --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --dosage 1 --keep_vcf filekeep_vcf

~/nextflow run gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_1vcf --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --dosage 1

~/nextflow run gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_multivcf  --listfile_vcf listvcf --dosage 1

~/nextflow run gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_multivcf_keepind  --listfile_vcf listvcf --dosage 1 --keep_vcf filekeep_vcf

~/nextflow run gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity  --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_multivcf  --listfile_vcf listvcf --dosage 1
```

