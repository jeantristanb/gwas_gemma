plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores

include {strmem} from './utils.nf'

process format_summarystat_gemmadosage{
 input :
     tuple val(pheno) ,path(filegemma)
 output :
     tuple val(pheno),path(filegemmaformat)
 //publishDir "${params.output_dir}/gemma/", mode:'copy'
 script :
     filegemmaformat=pheno+"_format.gemma"
     """
     format_gemmabimbam.py $filegemma $filegemmaformat
     """

}

process doGemmabimbam{
       label 'gemma'
       maxForks params.max_forks
       cpus params.gemma_num_cores
       time   params.big_time
       memory { strmem(params.gemma_mem_req) + 5.GB * (task.attempt -1) }
       errorStrategy { task.exitStatus in 137..144 ? 'retry' : 'terminate' }
       maxRetries 10
       input:
         tuple val(chro), path(rel_matrix),path(bimbam),path(bimbam_ind), path(annotation) ,path(data),path(rsfilelist), val(this_pheno), val(covariates), val(outdir)
       publishDir "${params.output_dir}/$outdir", mode:'copy'
       output:
         tuple val(our_pheno),path("${dir_gemma}/${out}.assoc.txt"), emit :resgemma
         path("${dir_gemma}/${out}.log.txt"), emit : log
       script:
          our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
          our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
          our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
          gemma_covariate    = "${our_pheno}.gemma_cov"
          phef               = "${our_pheno}_n.phe"
          covariate_option = (covariates=="") ? "" : " --cov_list ${covariates} " 
          bimbamhead=bimbam.baseName
          out                = "$our_pheno-$bimbamhead"
          covar_opt_gemma    =  (covariates=="") ? "" :  " -c $gemma_covariate " 
          dir_gemma          =  "gemma"
          annotation2        = (annotation=="") ? "" : " -a $annotation "

          """
          all_covariate.py --data  $data --bimbam_ind  $bimbam_ind $covariate_option --cov_out $gemma_covariate \
          --pheno $our_pheno2 --phe_out ${phef} --form_out 5
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          ${params.gemma_bin} -g $bimbam ${covar_opt_gemma}  -k $rel_matrix -lmm ${params.gemma_lmm}  -n 1 -p $phef -o $out -maf ${params.cut_maf} 
          mv output ${dir_gemma}
          """
}
process doGemma{
       label 'gemma'
       maxForks params.max_forks
       cpus params.gemma_num_cores
       memory { strmem(params.gemma_mem_req) + 5.GB * (task.attempt -1) }
       errorStrategy { task.exitStatus in 137..144 ? 'retry' : 'terminate' }
       maxRetries 10
       time   params.big_time
       input:
         tuple val(chro), path(rel),path(data), path(bed), path(bim), path(fam), path(rsfilelist), val(this_pheno), val(covariates),val(outdir)
       publishDir "${params.output_dir}/$outdir",  mode:'copy', pattern: '*.log'
       output:
         tuple val(our_pheno),path("${dir_gemma}/${out}.assoc.txt"), emit :resgemma
         path("${dir_gemma}/${out}.log.txt"), emit : log
       script:
          our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
          our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
          our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
          data_nomissing     = "pheno-"+our_pheno+".pheno"
          list_ind_nomissing = "lind-"+our_pheno+".lind"
          rel_matrix         = "newrel-"+our_pheno+".rel"
          base               =  bed.baseName
          inp_fam            =  base+".fam"
          newbase            =  base+"-"+our_pheno
          newfam             =  newbase+".fam"
          gemma_covariate    = "${newbase}.gemma_cov"
          phef               = "${newbase}_n.phe"
          covar_opt_gemma    =  (covariates) ?  " -c $gemma_covariate " : ""
          rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
          out                = "$base-$our_pheno-$chro"
          dir_gemma          =  "gemma"
          chroptionplk       =  (chro==-1) ? "" : "--chr $chro"
          covariate_option = (covariates) ?  " --cov_list ${covariates} " : ""
          """
          list_ind_nomissing.py --data $data --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
                                --lindout $list_ind_nomissing
          gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
          plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem} $chroptionplk
          all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                             --pheno $our_pheno2 --phe_out ${phef} --form_out 1
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm ${params.gemma_lmm}  -n 1 -p $phef -o $out -maf ${params.cut_maf} 
          mv output ${dir_gemma}
          rm $rel_matrix
          rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
          """
     }



process addNtoStatGemma{
           errorStrategy { task.exitStatus in 137..144 ? 'retry' : 'terminate' }
           memory { strmem(other_mem_req) + 5.GB * (task.attempt -1) }
          maxRetries 10
          input :
            tuple val(our_pheno), path(filestat), path(fileN)
          output :
            tuple val(our_pheno), path(newfilestat)
          script :
            newfilestat = "${our_pheno}_withN.gemma"
            """
            addn_statgwas.py --file_stat $filestat --file_freq $fileN --gwas_chr chr --gwas_ps ps --gwas_rs rs --out $newfilestat
            """
}

process doMergeGemma{
            input :
               tuple val(this_pheno),path(list_file)
            //publishDir "${params.output_dir}/gemma", mode:'copy'
            output :
                tuple val(our_pheno2), path("$out")
            script :
                our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
                our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
                out = "${our_pheno}.gemma"
                fnames = list_file.join(" ")
                file1  = list_file[0]
                """
                head -1 $file1 > $out
                cat $fnames | grep -v "p_wald" >> $out
                """
}

