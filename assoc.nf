#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *      Scott Hazelhurst
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2022
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths
nextflow.enable.dsl=2

  def newNamePheno(Pheno){
      SplP=Pheno.split(',')
      for (i = 0; i <SplP.size(); i++) {
         SplP[i]=(i+1)+"@@@"+SplP[i]
      }
      return(SplP)
  }




def helps = [ 'help' : 'help' ]

allowed_params = ["bfile","input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores","gemma_mem_req","gemma","linear","logistic","assoc","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "other_mem_req", "max_plink_cores", "pheno","big_time","thin", "gemma_mat_rel","print_pca", "listsnps_buildrelat","genetic_map_file", "rs_list","adjust","bootStorageSize","shared-storage-mount","mperm","sharedStorageMount","max-instances","maxInstances","boot-storage-size","sharedStorageMound","instance-type","instanceType","AMI", "gemma_multi",  "saige", 'gemma_loco', 'file_vcf', 'listfile_vcf', 'dosage', 'file_bimbam', 'file_bimbam_ind']

allowed_params_rel=["snps_exclude_rel", "snps_include_rel", "listsnps_buildrelat", "sample_snps_rel",  "thin_snp_rel", "cut_maf_rel"]
allowed_params+=allowed_params_rel


ParamBolt=["bolt_ld_scores_col", "bolt_ld_score_file","boltlmm", "bolt_covariates_type",  "bolt_use_missing_cov", "bolt_num_cores", "bolt_mem_req", "exclude_snps", "bolt_impute2filelist", "bolt_impute2fidiid", "bolt_otheropt","bolt_bin"]
allowed_params+=ParamBolt
ParamFast=["fastlmm","fastlmm_num_cores", "fastlmm_mem_req", "fastlmm_multi", "fastlmmc_bin", "covariates_type"]
allowed_params+=ParamFast
/*Gxe : */
GxE_params=['gemma_gxe', "plink_gxe", "gxe"]
allowed_params+=GxE_params
Saige_params=['pheno_bin', "list_vcf", "saige_num_cores", "saige_mem_req"]
allowed_params+=Saige_params

FastGWA_params=["fastgwa_mem_req", "fastgwa_num_cores", 'fastgwa', "gcta64_bin"]
allowed_params+=FastGWA_params

params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}



def params_help = new LinkedHashMap(helps)


filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/input"

params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
params.thin       = ""
params.covariates = ""
params.chrom      = ""
params.print_pca = 1
params.genetic_map_file = ""
params.list_vcf=""
params.vcf_field="DS"
params.min_scoreinfo=0.3
params.vcf_minmac=1
outfname = params.output_testing
params.cut_maf=0.01

params.listsnps_buildrelat = ""
params.snps_include_rel=""
params.snps_exclude_rel=""
params.sample_snps_rel=1
params.cut_maf_rel="0.01"
params.plink_indep_pairwise="100 20 0.1"
params.thin_snp_rel=""

/* Defines the path where any scripts to be executed can be found.
 */



/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 00

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests_all = ["assoc","fisher","model","cmh","linear","logistic","boltlmm", "fastlmm", "gemma", "gemma_gxe", 'saige']


params.assoc     = 0
params.fisher   = 0
params.cmh     =  0
params.model   =  0
params.linear   = 0
params.logistic = 0
params.gemma = 0
params.saige=0

params.gemma_multi=0
params.gemma_mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4
params.gemma_mat_rel = ""
params.gemma_num_cores = 8
params.gemma_loco = 0
params.file_bimbam = ""
params.listfile_bimbam = ""
params.listfile_vcf= ""

params.dosage = 0 
params.file_vcf=""

params.pheno = "_notgiven_"

//
params.saige_bin_fitmodel="/usr/local/bin/step1_fitNULLGLMM.R"
params.saige_bin_spatest="/usr/local/bin/step2_SPAtests.R"
params.saige_loco=1
params.saige_mem_req='10GB'
params.saige_num_cores=10

if (params.pheno == "_notgiven_") {
  println "No phenotype given -- set params.pheno";
  System.exit(-2);
}
  

/*JT Append initialisation variable*/
params.bolt_covariates_type = ""
params.bolt_ld_score_file= ""
params.bolt_ld_scores_col=""
params.boltlmm = 0
params.bolt_num_cores=8
params.bolt_mem_req="6GB"
params.bolt_use_missing_cov=0
params.exclude_snps=""
params.bolt_impute2filelist=""
params.bolt_impute2fidiid=""
params.bolt_otheropt=""
/*fastlmm param*/
params.bolt_bin ="bolt"
params.gemma_bin ="gemma"

/*gxe param : contains column of gxe*/
params.gemma_gxe=0
params.plink_gxe=0
params.max_plink_cores = 4
params.rs_list=""
params.gxe=""

/**/
params.fastgwa=0
params.fastgwa_mem_req="10G"
params.fastgwa_num_cores=5
params.grm_nbpart=100
params.grm_maf = 0.01
params.gcta64_bin = "gcta64"
params.fastgwa_type="--fastGWA-mlm-exact"
params.grm_cutoff =  0.05
params.covariates_type=""
params.gcta_grmfile=""
params.pheno_bin=0

params.genotype_field="GP"
params.score_imp="INFO"
params.qctoolsv2_bin="qctool"
params.bcftools_bin="bcftools"



params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '6GB' // how much plink needs for this
params.other_process_mem_req = '10G' // how much other processed need


plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores 
params.help = false
/*bfile definition*/
bfile=""
if(params.input_dir!="" && params.input_pat!=""){
 bfile=params.input_dir+"/"+params.input_pat
}else{
 if(params.bfile==""){
  println("bfile params or input_dir and output_dir not initialise")
 }else{
  bfile=params.bfile
 }
}
process plinkextractpos{
 cpus max_plink_cores
 input:
  tuple path(bed), path(bim), path(fam)
  path(listsnp)
  path(indkeep)
 output :
  tuple path("${bfilef}.bed"), path("${bfilef}.bim"), path("${bfilef}.fam")
 script :
   bfile=bed.baseName
   bfilef=bfile+'_relsub'
   """
   sed '1d' $indkeep | awk '{print \$1\"\\t\"\$2}' > indkeep
   plink -bfile $bfile  --extract range $listsnp --keep indkeep --keep-allele-order --make-bed -out $bfilef
   """
}

process plinkextractind{
 cpus max_plink_cores
 input :
  tuple path(bed), path (bim), path(fam)
  path(indkeep)
 output :
  tuple path("${outbed}.bed"), path("${outbed}.bim"), path("${outbed}.fam"), emit: filterind
 script :
  bfile=bed.baseName
  outbed=bfile+"_subind"
  """
  sed '1d' $indkeep | awk '{print \$1\"\\t\"\$2}' > indkeep
  plink -bfile $bfile --keep-allele-order --make-bed -out $outbed --threads $max_plink_cores --keep $indkeep
  """
}

process format_summarystat_gemmadosage{
 input :
     tuple val(pheno) ,path(filegemma)
 output :
     tuple val(pheno),path(filegemmaformat)
 //publishDir "${params.output_dir}/gemma/", overwrite:true, mode:'copy'
 script :
     filegemmaformat=pheno+"_format.gemma"
     """
     format_gemmabimbam.py $filegemma $filegemmaformat
     """

}

process getListeChro{
  input :
   tuple path(bed), path(bim), path(fam)
  output :
     stdout
  script:
     """
     cat $bim|awk '{print \$1}'|uniq|sort|uniq
     """
}

process formatvcfinbimbam{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     tuple val(chro), path(vcf)
  publishDir "${params.output_dir}/format/bimbam", overwrite:true, mode:'copy'
  output :
     tuple val(chro),path("${Ent}.bimbam"), path("${fileind}")
  script :
    headvcf=vcf.baseName
    Ent=(chro!=-1) ? "${headvcf}_${chro}" :  "$headvcf"
    chroparam=(chro!=-1) ?  " --regions $chro" : ""
    fileind=Ent+".ind"
    """
    zcat $vcf |head -10000|grep "#"|tail -1| awk '{for(Cmt=10;Cmt<=NF;Cmt++)print \$Cmt}' > $fileind
    bcftools index $vcf
    ${params.bcftools_bin} view -i '${params.score_imp}>${params.min_scoreinfo}' $chroparam $vcf |${params.qctoolsv2_bin} -g - -vcf-genotype-field ${params.genotype_field} -ofiletype bimbam_dosage -og ${Ent}.bimbam -filetype vcf
    """
}




process subsample_snps{
 cpus max_plink_cores
 input:
  tuple path(bed), path(bim), path(fam)
  path(snp_exclude_bed)
  path(snp_include_bed)
  path(indkeep)
 output:
  path("${outfile}.prune.in"), emit: subsample_snps_list
 script:
  bfile=bed.baseName
  outfile="indep_pairwise"
  rangeexclude=(params.snps_exclude_rel=="") ? "" : " --exclude range $snp_exclude_bed "
  rangeinclude=(params.snps_include_rel=="") ? "" : " --extract range $snp_include_bed "
  maf=(params.cut_maf_rel=="") ? "" : " --maf ${params.cut_maf_rel}"
  balisethin=(params.thin_snp_rel=="") ? "0" : "1"
  """
  sed '1d' $indkeep | awk '{print \$1\"\\t\"\$2}' > indkeep
  plink --bfile $bfile --threads $max_plink_cores --autosome $rangeexclude --indep-pairwise ${params.plink_indep_pairwise} --out $outfile $maf --keep indkeep --keep-allele-order $rangeinclude
  if [ "$balisethin" == "1" ]
  then
   cp ${outfile}.prune.in ${outfile}.prune.tmp.in
   shuf ${outfile}.prune.tmp.in | head -${params.thin_snp_rel} > ${outfile}.prune.in
  fi
  """
}
/*process insure that file is in chr pos pos rs*/
process checkposrsfile{
  input :
    path(pos)
    tuple path(bed), path(bim), path(fam) 
    val(out)
  output :
    path("$out")
  script :
     """
     check_filpos.py $pos $bim $out
     """
}

workflow getsnpbuilrelat{
 take :
   ch_plkfile
 main :
   checkposrsfile(channel.fromPath(params.listsnps_buildrelat, checkIfExists:true), ch_plkfile,'snpbuildrelat.bed')
 emit :
   pos_chr=checkposrsfile.out
}
workflow getsnpexcluderelat{
  take:
   ch_plkfile
  main:
  if(params.snps_exclude_rel!="")ch_snps_exclude_rel=checkposrsfile(channel.fromPath(params.snps_exclude_rel,checkIfExists:true),  ch_plkfile,'snpexcluderelat.bed').out
   else ch_snps_exclude_rel=channel.fromPath("${dummy_dir}/00")
  emit : 
   pos_chr=ch_snps_exclude_rel
}
workflow getsnpincluderelat{
  take :
   ch_plkfile
  main:
  if(params.snps_include_rel!="")ch_snps_include_rel=checkposrsfile(channel.fromPath(params.snps_include_rel,checkIfExists:true),  ch_plkfile,'snpincluderelat.bed').out
   else ch_snps_include_rel=channel.fromPath("${dummy_dir}/01")
  emit : 
   pos_chr=ch_snps_include_rel
}



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
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time params.big_time
       input:
         tuple path(bed), path(bim), path(fam)
       publishDir "${params.output_dir}/gemma/rel", overwrite:true, mode:'copy'
       output:
          tuple val(-1),path("output/${base}.*XX.txt") 
       script:
          base = bed.baseName
          famfile=base+".fam"
          """
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
          ${params.gemma_bin} -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 -km 2
          """
}

process getGemmaRelChro{
       label 'gemma'
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time params.big_time
       input:
         tuple path(bed), path(bim), path(fam), val(chro)
       publishDir "${params.output_dir}/gemma/rel", overwrite:true, mode:'copy'
       output:
          tuple val(chro), path("output/${newbase}.*XX.txt")
       script:
          base = bed.baseName
          newbase=base+"_${chro}"
          famfile=base+".fam"
          """
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
          plink -bfile $base --not-chr $chro --keep-allele-order --make-bed -out $newbase
          ${params.gemma_bin} -bfile $newbase  -gk ${params.gemma_relopt} -o $newbase -p pheno -n 3 -km 2 
          rm -rf $newbase*
          """
}

process GemmaBimbamRel{
       label 'gemma'
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time params.big_time
       input:
         tuple val(chro), path(bimbam), path(ind), path(listpos)
       publishDir "${params.output_dir}/gemma/rel", overwrite:true, mode:'copy'
       output:
          tuple val(chro),path("output/${base}.*XX.txt")
       script:
          tmp=bimbam.baseName
          base=(chro==-1) ? "${tmp}" : "${tmp}_${chro}"
          outposbimbam="newbimbam"
          """
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          cat $ind|awk '{print 0.2}' > pheno
          listpos_bimbam.py --bimbam $bimbam --filepos $listpos --out $outposbimbam --exclude_chr $chro
          ${params.gemma_bin} -g $outposbimbam -gk ${params.gemma_relopt} -o $base -p pheno -n 1 -km 1
          """
}


process doGemma{
       label 'gemma'
       maxForks params.max_forks
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time   params.big_time
       input:
         tuple val(chro), path(rel),path(covariates), path(bed), path(bim), path(fam), path(rsfilelist), val(this_pheno), val(outdir)
       publishDir "${params.output_dir}/$outdir", overwrite:true, mode:'copy', pattern: '*.log'
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
          covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
          rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
          out                = "$base-$our_pheno-$chro"
          dir_gemma          =  "gemma"
          chroptionplk       = 	(chro==-1) ? "" : "--chr $chro"
          covariate_option = (params.covariates) ?  " --cov_list ${params.covariates} " : "" 
          """
          list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
                                --lindout $list_ind_nomissing
          gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
          plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem} $chroptionplk
          all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                             --pheno $our_pheno2 --phe_out ${phef} --form_out 1
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf ${params.cut_maf}
          mv output ${dir_gemma}
          rm $rel_matrix
          rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
          """
     }
process doGemmabimbam{
       label 'gemma'
       maxForks params.max_forks
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time   params.big_time
       input:
         tuple val(chro), path(rel_matrix),path(bimbam),path(bimbam_ind),path(covariates),path(rsfilelist), val(this_pheno), val(outdir)
       publishDir "${params.output_dir}/$outdir", overwrite:true, mode:'copy', pattern: '*.log'
       output:
         tuple val(our_pheno),path("${dir_gemma}/${out}.assoc.txt"), emit :resgemma
         path("${dir_gemma}/${out}.log.txt"), emit : log
       script:
          our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
          our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
          our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
          gemma_covariate    = "${our_pheno}.gemma_cov"
          phef               = "${our_pheno}_n.phe"
          covariate_option = (params.covariates) ?  " --cov_list ${params.covariates} " : "" 
          bimbamhead=bimbam.baseName
          out                = "$our_pheno-$bimbamhead"
          covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
          dir_gemma          =  "gemma"

          """
          all_covariate.py --data  $covariates --bimbam_ind  $bimbam_ind $covariate_option --cov_out $gemma_covariate \
          --pheno $our_pheno2 --phe_out ${phef} --form_out 5
          export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
          ${params.gemma_bin} -g $bimbam ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf ${params.cut_maf}
          mv output ${dir_gemma}
          """
     }


     process doMergeGemma{
            input :
               tuple val(this_pheno),path(list_file)  
            //publishDir "${params.output_dir}/gemma", overwrite:true, mode:'copy'
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
    process splitbimbamchro{
        input :
          tuple val(chro), path(bimbam), path(bimbam_ind)
         output :
          tuple val(chro), path("$newfile"), path(bimbam_ind)
        script :
          newfile = bimbam.baseName.replaceAll(/.vcf$/,'')+"_" + chro+'.bimbam'
          """
          listpos_bimbam.py --bimbam $bimbam --include_chr $chro --out $newfile 
          """
    }
    process get_chrovcf{
       input:
          path(vcf)
       output :
           tuple env(chro), path(vcf), emit : chro_vcf
       script :
          """
          chro=`zcat $vcf|head -1000|grep -v "#"|awk '{print \$1}'|uniq`
          """
    }

process mergebimbamrel{
  input : 
    path(listbimam)
    path(listind)
    path(filepos) 
  output :
    tuple path(subbimbam), path("${listind[0]}")
  script :
    allbimbam=listbimam.join(',')
    subbimbam='allrelpos.bimbam'
    subbimbamnd='allrelpos.ind'
     """
     listpos_bimbam.py --listbimbam $allbimbam --filepos $filepos --out $subbimbam
     """
}

process getreport{
  input :
    tuple val(pheno), path(filegwas), val(outputdir)
  publishDir "${params.output_dir}/$outputdir/", overwrite:true, mode:'copy'
  output :
    file(filegwas) 
  script:
   """
   echo $filegwas
   """
}


workflow gwasgemma{
  take :
   ch_plkfile
   ch_plkfile_rel
   listchro
   filepheno
   filers
   listpheno
   bed_file_rel
 main :
   listchro_ch=listchro.flatMap{ list_str -> list_str.split() }
   if(params.dosage==1){
     if(params.file_bimbam!=""){
         bimbamfilei=channel.fromPath(params.file_bimbam,checkIfExists:true).combine(channel.fromPath(params.file_bimbam_ind, checkIfExists:true))
         bimbamfilerel=channel.fromPath(params.file_bimbam,checkIfExists:true).combine(channel.fromPath(params.file_bimbam_ind, checkIfExists:true))
         if(params.gemma_loco==1){
            splitbimbamchro(listchro_ch.combine(bimbamfilei))
            bimbamfile=splitbimbamchro.out
         }else{
             bimbamfile=bimbamfilei
         }
         
     }else if(params.file_vcf!=""){
         formatvcfinbimbam(channel.from(-1).combine(channel.fromPath(params.file_vcf, checkIfExists:true)))
         if(params.gemma_loco==1){
            splitbimbamchro(listchro_ch.combine(formatvcfinbimbam.out.flatMap{[it[1],it[2]]}.collect()))
           
            bimbamfile=splitbimbamchro.out
            bimbamfilerel=formatvcfinbimbam.out.flatMap{[it[1],it[2]]}.collect()
         }else{
         bimbamfile=formatvcfinbimbam.out.flatMap{[it[1],it[2]]}.collect()
         bimbamfilerel=formatvcfinbimbam.out.flatMap{[it[1],it[2]]}.collect()
         }
     }else if(params.listfile_bimbam!=""){
          chrobimbamfileI=channel.from(file(params.listfile_bimbam).readLines()).flatMap{it.split()[0]}
         //chrobimbamfileI.view()
          namebimbamfileI=channel.from(file(params.listfile_bimbam).readLines()).map{tuple(it.split()[0],file(it.split()[1]))}.combine(channel.fromPath(params.file_bimbam_ind, checkIfExists:true))
          bimbamfileI= namebimbamfileI//chrobimbamfileI.phase(namebimbamfileI)//.combine(channel.fromPath(params.file_bimbam_ind, checkIfExists:true))
          file_ch_bimbam=bimbamfileI.flatMap{it->it[1]}.collect()
          ind_ch_bimbam=channel.fromPath(params.file_bimbam_ind)
          mergebimbamrel(file_ch_bimbam, ind_ch_bimbam,bed_file_rel)
         if(params.gemma_loco==0)bimbamfile=bimbamfileI.flatMap{[it[1],it[2]].combinations()}
         else bimbamfile=bimbamfileI
         bimbamfilerel=mergebimbamrel.out
     }else if(params.listfile_vcf!=""){
         get_chrovcf(channel.fromPath(file(params.listfile_vcf, checkIfExists:true).readLines(), checkIfExists:true))
         formatvcfinbimbam(get_chrovcf.out.chro_vcf)
         file_ch_bimbam=formatvcfinbimbam.out.flatMap{it->it[1]}.collect()
         ind_ch_bimbam=formatvcfinbimbam.out.flatMap{it->it[2]}.collect()
         mergebimbamrel(file_ch_bimbam, ind_ch_bimbam,bed_file_rel)
         bimbamfilerel=mergebimbamrel.out
         if(params.gemma_loco==0)bimbamfile=formatvcfinbimbam.out.flatMap{[it[1],it[2]].combinations()}
         else bimbamfile=formatvcfinbimbam.out
     }else{
    println "No file gave for dosage, vcf imputation or bimbam file";
    System.exit(-2);

    }
     listpheno=Channel.from(listpheno)  
     if(params.gemma_loco==0){
         GemmaBimbamRel(channel.from("-1").combine(bimbamfilerel).combine(bed_file_rel))
//         tuple val(chro), path(rel_matrix),path(bimbam),path(bimbam_ind),path(covariates),path(rsfilelist), val(this_pheno), val(outdir)

         doGemmabimbam(GemmaBimbamRel.out.combine(bimbamfile).combine(filepheno).combine(filers).combine(listpheno).combine(channel.of('gemma/')))
         doMergeGemma(doGemmabimbam.out.resgemma.groupTuple())
     }else{
            GemmaBimbamRel(listchro_ch.combine(bimbamfilerel).combine(bed_file_rel))
            doGemmabimbam(GemmaBimbamRel.out.join(bimbamfile).combine(filepheno).combine(filers).combine(listpheno).combine(channel.of('gemma/chro')))
            doMergeGemma(doGemmabimbam.out.resgemma.groupTuple())
     }
    format_summarystat_gemmadosage(doMergeGemma.out)
    ressummstat=format_summarystat_gemmadosage.out
   }else{
    if(params.gemma_loco==0){
      getGemmaRelAll(ch_plkfile_rel)
      doGemma(getGemmaRelAll.out.combine(filepheno).combine(ch_plkfile).combine(filers).combine(listpheno).combine(channel.of('gemma/log/')))
      ressummstat=doGemma.out.resgemma
    }else{
     getGemmaRelChro(ch_plkfile_rel.combine(listchro_ch))
     doGemma(getGemmaRelChro.out.combine(filepheno).combine(ch_plkfile).combine(filers).combine(listpheno).combine(channel.of('gemma/chro/')))
     doMergeGemma(doGemma.out.resgemma.groupTuple())
     ressummstat=doMergeGemma.out
   }
 }
 getreport(ressummstat.combine(channel.of( 'gemma')))

}


workflow {
 /*bedfile*/
 bedfileI=Channel.fromPath("${bfile}.bed",checkIfExists:true).combine(Channel.fromPath("${bfile}.bim",checkIfExists:true)).combine(Channel.fromPath("${bfile}.fam",checkIfExists:true))
 if(params.rs_list=="")rsfile=Channel.fromPath("${dummy_dir}/06", checkIfExists:true)
 else rsfile=Channel.fromPath(params.rs_list, checkIfExists:true)
 phenofile=Channel.fromPath(params.data, checkIfExists:true)
 plinkextractind(bedfileI,phenofile)
 subsample_snp_rel(plinkextractind.out.filterind,phenofile)
 getListeChro(plinkextractind.out.filterind)
 listpheno = newNamePheno(params.pheno)
 gwasgemma(plinkextractind.out.filterind, subsample_snp_rel.out.plk_rel, getListeChro.out, phenofile, rsfile, listpheno, subsample_snp_rel.out.bed_pos_rel)

}

