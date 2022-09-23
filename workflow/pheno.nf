plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores

process add_pcs{
  cpus params.max_plink_cores
  input :
      tuple path(bed), path (bim), path(fam), path(phenofile)
  publishDir "${params.output_dir}/pheno_k/", mode:'copy'
  output :
    path(newfilepheno), emit :pheno
    env pcs, emit : pcs_covar
  script :
     plkf=bed.baseName
     head=phenofile.baseName
     newfilepheno=head+"_"+params.addpcs+".newpheno"

     """
     plink -bfile $plkf --pca ${params.addpcs} -keep $phenofile -out $head
     addpheno_pcs.r --data $phenofile --pcs $head".eigenvec" --out $newfilepheno
     pcs=`cat ${newfilepheno}.covar`
     """
}
 
workflow 
