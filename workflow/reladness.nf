plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores

include {getsnpexcluderelat} from './utils.nf'
include {getsnpincluderelat} from './utils.nf'
include {checkposrsfile} from './utils.nf'
include {subsample_snps} from './plink_utils.nf'
include {plinkextractpos} from './plink_utils.nf'
include {strmem} from './utils.nf'



workflow subsample_snp_rel{
 take :
  ch_plkfile
  ch_pheno
 main :
  if(params.listsnps_buildrelat!=""){
    getsnpbuilrelat(ch_plkfile)
    snpfilers=getsnpbuilrelat.out.posfile
  }else{
   getsnpexcluderelat(ch_plkfile)
   getsnpincluderelat(ch_plkfile)
   subsample_snps(ch_plkfile, getsnpexcluderelat.out.pos_chr,getsnpincluderelat.out.pos_chr,ch_pheno)
   checkposrsfile(subsample_snps.out.subsample_snps_list, ch_plkfile,'list_posrs_rel')
   snpfilers=checkposrsfile.out
  }
 plinkextractpos(ch_plkfile, snpfilers, ch_pheno)
 emit:
   bed_pos_rel=snpfilers
   plk_rel=plinkextractpos.out
}


process getGemmaRelAll {
       label 'gemma'
       cpus params.gemma_num_cores_rel
       maxForks params.max_forks
       memory { strmem(params.gemma_mem_req_rel) + 5.GB * (task.attempt -1) }
       errorStrategy { task.exitStatus in 137..144 ? 'retry' : 'terminate' }
       maxRetries 10
       time params.big_time
       input:
         tuple path(bed), path(bim), path(fam)
       publishDir "${params.output_dir}/gemma/rel", mode:'copy'
       output:
          tuple val(-1),path("output/${base}.*XX.txt"), emit : rel
          path("output/*.log.txt"), emit : log
       script:
          base = bed.baseName
          famfile=base+".fam"
          """
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores_rel}
          cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
          ${params.gemma_bin} -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 -km 2
          """
}

process getGemmaRelChro{
       label 'gemma'
       cpus params.gemma_num_cores_rel
       memory { strmem(params.gemma_mem_req_rel) + 5.GB * (task.attempt -1) }
       errorStrategy { task.exitStatus in 137..144 ? 'retry' : 'terminate' }
       maxForks params.max_forks
       maxRetries 10
       time params.big_time
       input:
         tuple path(bed), path(bim), path(fam), val(chro)
       publishDir "${params.output_dir}/gemma/rel", mode:'copy'
       output:
          tuple val(chro), path("output/${newbase}.*XX.txt"), emit :rel
          path("output/*.log.txt"), emit : log
       script:
          base = bed.baseName
          newbase=base+"_${chro}"
          famfile=base+".fam"
          """
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores_rel}
          cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
          plink -bfile $base --not-chr $chro --keep-allele-order --make-bed -out $newbase
          ${params.gemma_bin} -bfile $newbase  -gk ${params.gemma_relopt} -o $newbase -p pheno -n 3 -km 2
          rm -rf $newbase*
          """
}


process GemmaBimbamRel{
       label 'gemma'
       cpus params.gemma_num_cores_rel
       time params.big_time
       memory { strmem(params.gemma_mem_req_rel) + 5.GB * (task.attempt -1) }
       maxForks params.max_forks
       errorStrategy 'retry'
       maxRetries 10
       input:
         tuple val(chro), path(bimbam), path(ind),path(annotation),path(listpos)
       publishDir "${params.output_dir}/gemma/rel", mode:'copy'
       output:
          tuple val(chro),path("output/${base}.*XX.txt"), emit :rel
          path("output/*.log.txt"), emit : log
       script:
          tmp=bimbam.baseName
          base=(chro==-1) ? "${tmp}" : "${tmp}_${chro}"
          outposbimbam="posbimbam_"+chro
          annot="sub_annot_"+chro+".txt"
          """
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores_rel}
          cat $ind|awk '{print 0.2}' > pheno
          listpos_bimbam.py --bimbam $bimbam --filepos $listpos --out $outposbimbam --exclude_chr $chro --annotation $annot
          ${params.gemma_bin} -g $outposbimbam -gk ${params.gemma_relopt} -o $base -p pheno -n 1 -km 1 -a $annot
          """
}



