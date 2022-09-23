#nextflow pull h3abionet/h3agwas
## first test : gemma with plink no loco
ls ~/Travail/git/h3agwas_buildexample/data2/vcf/*.vcf.gz > listvcf
FilePheno=exampledata2b.pheno
#~/Travail/git/h3agwas_buildexample/data2//exampledata2.pheno
#awk '{$5=$4;print $0}' ~/Travail/git/h3agwas_buildexample/data2//exampledata2.pheno > exampledata2b.pheno

#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_plk_vcf  --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --reffasta ~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz 
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_plk_multivcf  --listfile_vcf listvcf --reffasta ~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz 

#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_1vcf --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --dosage 1
#zcat ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz | head -1000 |grep "#"  | tail -1 |awk '{for(Cmt=40;Cmt<NF;Cmt++)print $Cmt}' > filekeep_vcf
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_1vcf_indkeep --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --dosage 1 --keep_vcf filekeep_vcf
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_1vcf --file_vcf ~/Travail/git/h3agwas_buildexample/data2/exampledata2.result.imp.vcf.gz --dosage 1
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_multivcf  --listfile_vcf listvcf --dosage 1
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_multivcf_keepind  --listfile_vcf listvcf --dosage 1 --keep_vcf filekeep_vcf

#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_multivcf  --listfile_vcf listvcf --dosage 1

filebimbam=pheno_q1_loco_bimbam_1vcf/format/bimbam/exampledata2.result.imp.vcf.bimbam
filebimbamind=pheno_q1_loco_bimbam_1vcf/format/bimbam/exampledata2.result.imp.vcf.ind
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_1bimbam --file_bimbam $filebimbam --dosage 1 --file_bimbam_ind $filebimbamind
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_1bimbam --file_bimbam $filebimbam --dosage 1 --file_bimbam_ind $filebimbamind
#ls pheno_q1_loco_bimbam_multivcf/format/bimbam/*.bimbam > listbimbam
listfilebimbam=listbimbam
#rm -f  $listfilebimbam
#for chro in `seq 1 22`
#do
#filebimbam=$PWD/pheno_q1_loco_bimbam_multivcf/format/bimbam/"tmp_"$chro".pbwt_reference_impute.vcf_"$chro".bimbam"
#echo "$chro $filebimbam" >> $listfilebimbam
#done

#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_multibimbam --listfile_bimbam $listfilebimbam --dosage 1 --file_bimbam_ind $filebimbamind
#~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf --input_dir ~/Travail/git/h3agwas_buildexample/data2//qc/qc/ --input_pat exampledata2_qc --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 0   --output_dir pheno_q1_noloco_bimbam_multibimbam --listfile_bimbam $listfilebimbam --dosage 1 --file_bimbam_ind $filebimbamind


~/nextflow run /home/jeantristan/Travail/GWAS/GWAS_CKD/Version22/gwas_ckdawigen/assoc.nf  --data $FilePheno --pheno pheno_1,pheno_2  --output pheno1 --gemma 1 -profile slurmSingularity --assoc 1 --sample_snps_rel 1 --linear 1 -resume --gemma_loco 1   --output_dir pheno_q1_loco_bimbam_multivcf  --listfile_vcf listvcf --dosage 1 --reffasta ~/Travail/git/h3agwas_buildexample/data2/formatimput/hg19.fa.gz 
