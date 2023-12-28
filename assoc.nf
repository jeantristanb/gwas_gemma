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

allowed_params = ["bfile","input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores", "gemma_num_cores_rel","gemma_mem_req","gemma","linear","logistic","assoc","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "other_mem_req", "max_plink_cores", "pheno","big_time","thin", "gemma_mat_rel","print_pca", "listsnps_buildrelat","genetic_map_file", "rs_list","adjust","bootStorageSize","shared-storage-mount","sharedStorageMount","max-instances","maxInstances","boot-storage-size","sharedStorageMound","instance-type","instanceType","AMI", "gemma_multi",  "saige", 'gemma_loco', 'file_vcf', 'listfile_vcf', 'dosage', 'file_bimbam', 'file_bimbam_ind','keep_vcf', 'vcftools_bin', "gemma_lmm"]

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
params.input_dir  = ""
params.output_dir = "${params.work_dir}/output"
params.output_pat= "output"
params.bfile= ""
params.output_testing = "cleaned"
params.thin       = ""
params.covariates = ""
params.chrom      = ""
params.print_pca = 1
params.genetic_map_file = ""
params.list_vcf=""
params.min_scoreinfo=0.3
outfname = params.output_testing
params.cut_maf=0.01
params.keep_vcf=''
params.statfreq_vcf="%AN %AC"

params.listsnps_buildrelat = ""
params.snps_include_rel=""
params.snps_exclude_rel=""
params.sample_snps_rel=1
params.cut_maf_rel="0.01"
params.plink_indep_pairwise="100 20 0.1"
params.thin_snp_rel=""
params.vcfftools_bin='vcftools'

params.pheno_tr_fct="" 
params.phenores_tr_fct="" 
params.pheno_residuals=1
params.addpcs = 0

params.gemma_lmm=4

/* Defines the path where any scripts to be executed can be found.
 */



/* Do permutation testing -- 0 for none, otherwise give number */

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests_all = ["assoc","fisher","model","cmh","linear","logistic","boltlmm", "fastlmm", "gemma", "gemma_gxe", 'saige']


params.gemma = 0

params.gemma_multi=0
params.gemma_mem_req = "6GB"
params.gemma_mem_req_rel = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4
params.gemma_mat_rel = ""
params.gemma_num_cores = 8
params.gemma_num_cores_rel = 8
params.gemma_loco = 0
params.file_bimbam = ""
params.listfile_bimbam = ""
params.listfile_bimbam_annot = ""
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
if(params.gemma_mem_req_rel==""){
params.gemma_mem_req_rel = params.gemma_mem_req
}


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


/*params for format*/
params.file_listvcf=""


params.genetic_maps=""
params.do_stat=true
params.unzip_zip=0
params.unzip_password=""
params.reffasta=""
params.other_process_mem_req="10GB"



params.input_pat  = ''

params.sexinfo_available = "false"


params.plink_mem_req = '6GB' // how much plink needs for this
params.other_process_mem_req = '10G' // how much other processed need


plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores 
params.help = false
/*bfile definition*/

include {format_vcfinplk} from './workflow/convert_file.nf'
include {wf_prepare_pheno} from './workflow/pheno.nf'
include {splitbimbamchro} from './workflow/utils.nf'
include {formatvcfinbimbam} from './workflow/convert_file.nf'
include {formatvcfinbimbam_ind} from './workflow/convert_file.nf'
//include {mergebimbamrel} from './workflow/utils.nf'
include {mergebimbamrel_speed as mergebimbamrel} from './workflow/utils.nf'
include {strmem} from './workflow/utils.nf'
include {get_chrovcf} from './workflow/vcf.nf'
include {GemmaBimbamRel} from './workflow/reladness.nf'
include {getGemmaRelAll} from './workflow/reladness.nf'
include {getGemmaRelChro} from './workflow/reladness.nf'
include {doGemmabimbam} from './workflow/gemma.nf'
include {doMergeGemma} from './workflow/gemma.nf'
include {doGemma} from './workflow/gemma.nf'
include {addNtoStatGemma} from './workflow/gemma.nf'
include {computeN_plink} from './workflow/plink_utils.nf'
include {plinkextractind} from './workflow/plink_utils.nf'
include {subsample_snp_rel} from './workflow/reladness.nf'
include {getListeChro} from './workflow/utils.nf'
include {cleanvcfwf} from './workflow/vcf.nf'
include {format_summarystat_gemmadosage} from './workflow/gemma.nf'









process getreport{
  input :
    tuple val(pheno), path(filegwas), val(outputdir)
  publishDir "${params.output_dir}/$outputdir/", mode:'copy'
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
	   filevcf
	   listfilevcf
	   covariates
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
		 //formatvcfinbimbam(channel.from(-1).combine(filevcf))
		 formatvcfinbimbam_ind(channel.from(-1).combine(filevcf).combine(filepheno))
		 if(params.gemma_loco==1){
		    splitbimbamchro(listchro_ch.combine(formatvcfinbimbam_ind.out.bimbam.flatMap{[it[1],it[2],it[3]]}.collect()))
		    bimbamfile=splitbimbamchro.out
		    bimbamfilerel=formatvcfinbimbam_ind.out.bimbam.flatMap{[it[1],it[2], it[3]]}.collect()
		 }else{
		 bimbamfile=formatvcfinbimbam_ind.out.bimbam.flatMap{[it[1],it[2], it[3]]}.collect()
		 bimbamfilerel=formatvcfinbimbam_ind.out.bimbam.flatMap{[it[1],it[2], it[3]]}.collect()
		 }
	     }else if(params.listfile_bimbam!=""){
                  println("used file bimbam give by user")
                  if(params.listfile_bimbam_annot!=""){
		    bimbam_annot=channel.from(file(params.listfile_bimbam_annot).readLines()).map{tuple(it.split()[0],file(it.split()[1]))}
                  }else{
		    bimbam_annot=channel.from(file(params.listfile_bimbam).readLines()).flatMap{it.split()[0]}.combine(channel.fromPath("${dummy_dir}/07"))
                  }
		  chrobimbamfileI=channel.from(file(params.listfile_bimbam).readLines()).flatMap{it.split()[0]}
		  namebimbamfileI=channel.from(file(params.listfile_bimbam).readLines()).map{tuple(it.split()[0],file(it.split()[1]))}.combine(channel.fromPath(params.file_bimbam_ind, checkIfExists:true)).join(bimbam_annot)
		  bimbamfileI= namebimbamfileI//chrobimbamfileI.phase(namebimbamfileI)//.combine(channel.fromPath(params.file_bimbam_ind, checkIfExists:true))
		   file_ch_bimbam=bimbamfileI.flatMap{it->it[1]}.collect()
		   ind_ch_bimbam=channel.fromPath(params.file_bimbam_ind)
		   mergebimbamrel(bimbamfileI.flatMap{it->it[1]}, ind_ch_bimbam,bed_file_rel)
		 if(params.gemma_loco==0)bimbamfile=bimbamfileI.flatMap{[it[1],it[2], it[3]].combinations()}
		 else bimbamfile=bimbamfileI
		 bimbamfilerel=mergebimbamrel.out
	     }else if(params.listfile_vcf!=""){
		 get_chrovcf(listfilevcf)
		 //formatvcfinbimbam(get_chrovcf.out.chro_vcf)
		 formatvcfinbimbam_ind(get_chrovcf.out.chro_vcf.combine(filepheno))
		 file_ch_bimbam=formatvcfinbimbam_ind.out.bimbam.flatMap{it->it[1]}.collect()
		 ind_ch_bimbam=formatvcfinbimbam_ind.out.bimbam.flatMap{it->it[2]}.collect()
		 mergebimbamrel(formatvcfinbimbam_ind.out.bimbam.flatMap{it->it[1]}, ind_ch_bimbam,bed_file_rel)
		 bimbamfilerel=mergebimbamrel.out
		 if(params.gemma_loco==0)bimbamfile=formatvcfinbimbam_ind.out.flatMap{[it[1],it[2], it[3]].combinations()}
		 else bimbamfile=formatvcfinbimbam_ind.out.bimbam
	     }else{
	    println "No file gave for dosage, vcf imputation or bimbam file";
	    System.exit(-2);

	    }
	     if(params.gemma_loco==0){
		 GemmaBimbamRel(channel.from("-1").combine(bimbamfilerel).combine(bed_file_rel))
		 doGemmabimbam(GemmaBimbamRel.out.rel.combine(bimbamfile).combine(filepheno).combine(filers).combine(listpheno).combine(covariates).combine(channel.of('gemma/')))
		 doMergeGemma(doGemmabimbam.out.resgemma.groupTuple())
	     }else{
		    GemmaBimbamRel(listchro_ch.combine(bimbamfilerel).combine(bed_file_rel))
		    doGemmabimbam(GemmaBimbamRel.out.rel.join(bimbamfile).combine(filepheno).combine(filers).combine(listpheno).combine(covariates).combine(channel.of('gemma/chro')))
		    doMergeGemma(doGemmabimbam.out.resgemma.groupTuple())
	     }
	    format_summarystat_gemmadosage(doMergeGemma.out)
	    ressummstat=format_summarystat_gemmadosage.out
	   }else{
	    if(params.gemma_loco==0){
	      getGemmaRelAll(ch_plkfile_rel)
	      doGemma(getGemmaRelAll.out.rel.combine(filepheno).combine(ch_plkfile).combine(filers).combine(listpheno).combine(covariates).combine(channel.of('gemma/log/')))
	      ressummstat=doGemma.out.resgemma
	    }else{
	     getGemmaRelChro(ch_plkfile_rel.combine(listchro_ch))
	     doGemma(getGemmaRelChro.out.rel.combine(filepheno).combine(ch_plkfile).combine(filers).combine(listpheno).combine(covariates).combine(channel.of('gemma/chro/')))
	     doMergeGemma(doGemma.out.resgemma.groupTuple())
	     ressummstat=doMergeGemma.out
	   }
	 }
	 computeN_plink(filepheno.combine(ch_plkfile).combine(listpheno).combine(covariates))
	 addNtoStatGemma(ressummstat.join(computeN_plink.out))
	 getreport(addNtoStatGemma.out.combine(channel.of( 'gemma')))

}





workflow {
      bfile=""
      if(params.input_dir!="" && params.input_pat!="") bfile=params.input_dir+"/"+params.input_pat 
      else{
         if(params.bfile=="" && (params.file_vcf!="" || params.listfile_vcf!="")){
          println("bfile params or input_dir and output_dir not initialise")
          bfile=params.bfile
         } else bfile=params.bfile
        }
        println bfile
       //sys,exit()
	 /*bedfile*/
	 if(bfile!=""){
	   bedfileI=Channel.fromPath("${bfile}.bed",checkIfExists:true).combine(Channel.fromPath("${bfile}.bim",checkIfExists:true)).combine(Channel.fromPath("${bfile}.fam",checkIfExists:true))
	 }else{
	  println "no input plink format, used vcf file(s) to format in plink"
	  if(params.reffasta=="" || params.reffasta==true){
	   println("to format vcf in plink need a fasta file : args --reffasta null")
	      exit 1
	  }
          if(params.file_vcf=="" && params.listfile_vcf==""){
	      println("error : if no file vcf allowed (--file_vcf) or (--listfile_vcf), plink file must initalise (--input_dir / input_pat or --bfile)")
	      exit 10
          }
	  format_vcfinplk()
	  bedfileI=format_vcfinplk.out.plk
	 }
	 if(params.rs_list=="")rsfile=Channel.fromPath("${dummy_dir}/06", checkIfExists:true)
	 else rsfile=Channel.fromPath(params.rs_list, checkIfExists:true)
	 phenofile=Channel.fromPath(params.data, checkIfExists:true)
	 plinkextractind(bedfileI,phenofile)
	 subsample_snp_rel(plinkextractind.out.filterind,phenofile)
	 getListeChro(plinkextractind.out.filterind)
         wf_prepare_pheno(phenofile, subsample_snp_rel.out.plk_rel)
	 //listpheno = newNamePheno(params.pheno)
         listpheno=wf_prepare_pheno.out.pheno.flatMap{it->it.split(',')}
	 cleanvcfwf()
	 gwasgemma(plinkextractind.out.filterind, subsample_snp_rel.out.plk_rel, getListeChro.out,  wf_prepare_pheno.out.data, rsfile,listpheno, subsample_snp_rel.out.bed_pos_rel, cleanvcfwf.out.filevcf, cleanvcfwf.out.listfilevcf, wf_prepare_pheno.out.covar)

	}

