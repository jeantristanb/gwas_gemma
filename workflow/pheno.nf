plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores

process add_pcs{
  label 'R' 
  cpus params.max_plink_cores
  input :
      tuple path(bed), path (bim), path(fam), path(phenofile)
  publishDir "${params.output_dir}/pheno/pca/", mode:'copy'
  output :
    path(newfilepheno), emit :data
    env pcs, emit : pcs_covar
  script :
     plkf=bed.baseName
     head=phenofile.baseName
     newfilepheno=head+"_"+params.addpcs+".newpheno"
     covar2=(params.covariates=="")? "" : " --covar $params.covariates " 

     """
     plink -bfile $plkf --pca ${params.addpcs} -keep $phenofile -out $head
     addpheno_pcs.r --data $phenofile --pcs $head".eigenvec" --out $newfilepheno $covar2
     pcs=`cat ${newfilepheno}.covar`
     """
}

process format_pheno {
  label 'R'
  input :
       path(phenofile)
       tuple path(bed), path (bim), path(fam)
       val(pheno)
       val(covar)
  publishDir "${params.output_dir}/pheno/", mode:'copy'
  output :
       path("${newdatafile}.pheno"), emit :data
       env newpheno, emit : pheno
       val(newcovar), emit : covar
       path("${newdatafile}*")
  script :
      newdatafile=phenofile.baseName+'_format'
      covar2=(covar=="")? "" : " --covar $covar " 
      phenores_tr_fct= (params.phenores_tr_fct=="") ? "none"  : "$params.phenores_tr_fct"
      pheno_tr_fct= (params.pheno_tr_fct=="") ? "none"  : "$params.pheno_tr_fct"
      newcovar=(params.pheno_residuals==1) ? "" : " $covar" 
  """
  format_pheno.r --pheno ${pheno} --data $phenofile $covar2 --out $newdatafile --fam $fam --transform_i ${pheno_tr_fct} --transform_r ${phenores_tr_fct} --residuals ${params.pheno_residuals}
  newpheno=`cat pheno.txt`
  """
}
 
workflow wf_prepare_pheno{
   take : 
        phenofile 
        plink
   main :
    pheno=channel.from(params.pheno)
    covar=channel.from(params.covariates)
    if(params.addpcs>=1){
      add_pcs(plink.combine(phenofile))  
      new_data=add_pcs.out.data
      newcovar=add_pcs.out.pcs_covar
      //if(params.covar!='')newcovar=params.covar+","+newcovar
      data_ch=add_pcs.out.data
    }else {
    newcovar=params.covariates 
    data_ch=phenofile
   }
  format_pheno(data_ch, plink,pheno,newcovar) 
  emit :
     data=format_pheno.out.data 
     pheno=format_pheno.out.pheno
     covar=format_pheno.out.covar
}
