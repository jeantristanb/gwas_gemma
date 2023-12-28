
filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/input"

plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores


def strmem(val){
 return val as nextflow.util.MemoryUnit
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

/*check snp*/
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

  if(params.snps_include_rel!=""){
     checkposrsfile(channel.fromPath(params.snps_include_rel,checkIfExists:true),  ch_plkfile,'snpincluderelat.bed')
     ch_snps_include_rel=checkposrsfile.out
   }else ch_snps_include_rel=channel.fromPath("${dummy_dir}/01")
  emit :
   pos_chr=ch_snps_include_rel
}

process splitbimbamchro{
        input :
          tuple val(chro), path(bimbam), path(bimbam_ind), path(annotation)
         output :
          tuple val(chro), path("$newfile"), path(bimbam_ind),  path(annotation)
        script :
          newfile = bimbam.baseName.replaceAll(/.vcf$/,'')+"_" + chro+'.bimbam.gz'
          annotation= bimbam.baseName.replaceAll(/.vcf$/,'')+"_" + chro+'.annotation'
          """
          listpos_bimbam.py --bimbam $bimbam --include_chr $chro --out $newfile --annotation $annotation
          """
}

process mergebimbamrel{
  input :
    path(listbimam)
    path(listind)
    path(filepos)
  output :
    tuple path(subbimbam), path("listind.bimbam.out"), path(annotation)
  script :
    allbimbam=listbimam.join(',')
    subbimbam='allrelpos.bimbam.gz'
    subbimbamnd='allrelpos.ind'
    annotation="annot_rel.txt"
    """
     cp ${listind[0]} listind.bimbam.out
     listpos_bimbam.py --listbimbam $allbimbam --filepos $filepos --out $subbimbam --annotation $annotation
    """
}

process extract_bimbam_chr{
  input :
    tuple path(bimbam), path(listind), path(filepos)
  output :
    path(subbimbam), emit : bimbam
    path(subbimbamnd), emit : ind
    path(annotation), emit : annot
  script :
    subbimbam=bimbam+"_sub.bimbam.gz"
    subbimbamnd=bimbam+"_sub.ind"
    annotation=bimbam+"_sub.annotation"
    """
     cp ${listind[0]} $subbimbamnd
     chro=`zcat $bimbam|head -1|awk '{print \$1}'|awk -F":" '{print \$2}'`
     listpos_bimbam.py --listbimbam $bimbam --filepos $filepos --out $subbimbam --annotation $annotation --include_chr \$chro
    """
}

process merge_bimbam {
  input :
     path(bimbam)
     path(annot)
     path(listind)
  output :
    tuple path(subbimbam), path(subbimbamind), path(annotation)
  script :
    subbimbam='allrelpos.bimbam.gz'
    subbimbamind="listind.bimbam.out"
    annotation="annot_rel.txt"
    """
    cat *_sub.bimbam.gz > allrelpos.bimbam.gz
    cat *.annotation > $annotation
    cp ${listind[0]} $subbimbamind
    """
}

workflow mergebimbamrel_speed {
 take :
    listbimbam
    listind
    filepos
 main :
   extract_bimbam_chr(listbimbam.combine(listind).combine(filepos))
   merge_bimbam(extract_bimbam_chr.out.bimbam.collect(), extract_bimbam_chr.out.annot.collect(),extract_bimbam_chr.out.ind.collect())
 emit :
   merge_bimbam.out   

}
